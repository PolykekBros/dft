import os
import re
import numpy as np
from ase import Atoms, units
from ase.io import write
from ase.build import surface
from ase.constraints import FixAtoms

# --- Configuration ---
BULK_OUT = "../../volume/NiO_struct_opt_dft_U/NiO.vc_relax_dft_U.out"
PSEUDO_DIR = "../../../pot"
MILLER_INDICES = (1, 0, 0)
LAYERS_LIST = [4, 6, 8]
VACUUM_A = 14.0
FIX_BOTTOM_FRACTION = 0.5
SLAB_KPTS = (4, 4, 1)


# --- Custom QE Output Parser (Noncollinear/SOC Safe) ---
def load_qe_out(filepath):
    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()

    total_energy_ry = None
    cell = None
    species_labels = []
    positions = []
    pos_unit = "crystal"

    for i, line in enumerate(lines):
        if line.strip().startswith("!    total energy"):
            total_energy_ry = float(line.split("=")[1].split()[0])

        if line.strip().startswith("CELL_PARAMETERS"):
            match = re.search(r"alat\s*=\s*([0-9.+-eEdD]+)", line)
            if match:
                scale = float(match.group(1)) * units.Bohr
            elif "bohr" in line.lower():
                scale = units.Bohr
            else:
                scale = 1.0
            cell = (
                np.array(
                    [list(map(float, lines[i + k].split()[:3])) for k in (1, 2, 3)]
                )
                * scale
            )

        if line.strip().startswith("ATOMIC_POSITIONS"):
            if "crystal" in line.lower():
                pos_unit = "crystal"
            elif "bohr" in line.lower():
                pos_unit = "bohr"
            else:
                pos_unit = "angstrom"

            species_labels = []
            positions = []
            j = i + 1
            while j < len(lines) and lines[j].strip():
                parts = lines[j].split()
                if len(parts) < 4 or parts[0].startswith("End"):
                    break
                species_labels.append(parts[0])
                positions.append(list(map(float, parts[1:4])))
                j += 1

    symbols = [re.sub(r"\d+$", "", lbl) for lbl in species_labels]
    tags = [1 if lbl == "Ni1" else (2 if lbl == "Ni2" else 0) for lbl in species_labels]

    if pos_unit == "crystal":
        atoms = Atoms(symbols=symbols, scaled_positions=positions, cell=cell, pbc=True)
    else:
        scale = units.Bohr if pos_unit == "bohr" else 1.0
        atoms = Atoms(
            symbols=symbols, positions=np.array(positions) * scale, cell=cell, pbc=True
        )

    atoms.set_tags(tags)
    e_total_ev = total_energy_ry * units.Ry
    return atoms, e_total_ev


# --- QE 7.4.1 Noncollinear + SOC + DFT+U Input Writer ---
def write_qe_slab_input(filename, slab, fixed_mask, prefix, kpts, pseudo_dir):
    species_labels = []
    for sym, tag in zip(slab.get_chemical_symbols(), slab.get_tags()):
        if sym == "Ni":
            species_labels.append("Ni1" if tag == 1 else "Ni2")
        else:
            species_labels.append("O")

    lines = [
        "&CONTROL",
        "    calculation   = 'relax'",
        f"    prefix        = '{prefix}'",
        "    restart_mode  = 'from_scratch'",
        "    outdir        = './out'",
        f"    pseudo_dir    = '{pseudo_dir}'",
        "    tprnfor       = .true.",
        "    etot_conv_thr = 1e-5",
        "    forc_conv_thr = 1e-4",
        "/",
        "&SYSTEM",
        "    ibrav = 0",
        f"    nat   = {len(slab)}",
        "    ntyp  = 3",
        "    noncolin = .true.",
        "    lspinorb = .true.",
        "    starting_magnetization(1) = 0.5",
        "    starting_magnetization(2) = 0.5",
        "    angle1(1) = 0.0",
        "    angle1(2) = 180.0",
        "    angle2(1) = 0.0",
        "    angle2(2) = 0.0",
        "    ecutwfc = 65",
        "    ecutrho = 260",
        "    occupations = 'smearing'",
        "    smearing    = 'cold'",
        "    degauss     = 0.01",
        "    tefield  = .true.",
        "    dipfield = .true.",
        "    edir     = 3",
        "    emaxpos  = 0.95",
        "    eopreg   = 0.05",
        "    eamp     = 0.0",
        "/",
        "&ELECTRONS",
        "    mixing_beta      = 0.1",
        "    electron_maxstep = 300",
        "    mixing_ndim      = 14",
        "    conv_thr         = 1d-9",
        "/",
        "&IONS",
        "    ion_dynamics = 'bfgs'",
        "/",
        "ATOMIC_SPECIES",
        " Ni1  58.6934  Ni.upf",
        " Ni2  58.6934  Ni.upf",
        " O    15.999   O.upf",
        "",
        "CELL_PARAMETERS (angstrom)",
    ]

    for vec in slab.cell:
        lines.append(f"  {vec[0]:16.10f}  {vec[1]:16.10f}  {vec[2]:16.10f}")

    lines.append("")
    lines.append("ATOMIC_POSITIONS (angstrom)")
    for label, pos, is_fixed in zip(species_labels, slab.positions, fixed_mask):
        flag = "0 0 0" if is_fixed else "1 1 1"
        lines.append(
            f"  {label:<4s}  {pos[0]:14.8f}  {pos[1]:14.8f}  {pos[2]:14.8f}  {flag}"
        )

    lines.extend(
        [
            "",
            "K_POINTS automatic",
            f"{kpts[0]} {kpts[1]} {kpts[2]} 0 0 0",
            "",
            "HUBBARD (ortho-atomic)",
            "U Ni1-3d 5",
            "U Ni2-3d 5",
            "U O-2p 6",
            "",
        ]
    )

    with open(filename, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


# --- Post-Relaxation Surface Energy Verification ---
def verify_slab_convergence(bulk_out_path, slab_out_paths):
    bulk_ref, e_bulk_tot = load_qe_out(bulk_out_path)
    e_bulk_atom = e_bulk_tot / len(bulk_ref)

    for path in slab_out_paths:
        slab_rel, e_slab = load_qe_out(path)
        area = np.linalg.norm(np.cross(slab_rel.cell[0], slab_rel.cell[1]))
        gamma_ev_a2 = (e_slab - len(slab_rel) * e_bulk_atom) / (2.0 * area)
        gamma_j_m2 = gamma_ev_a2 * units._e / 1e-20
        print(f"{path}: gamma = {gamma_ev_a2:.4f} eV/A^2 ({gamma_j_m2:.4f} J/m^2)")


# --- Main Slab Generation Workflow ---
if __name__ == "__main__":
    bulk, e_bulk_total = load_qe_out(BULK_OUT)
    e_bulk_per_atom = e_bulk_total / len(bulk)
    print(
        f"Bulk formula: {bulk.get_chemical_formula()} | E_bulk: {e_bulk_per_atom:.6f} eV/atom"
    )

    for n_layers in LAYERS_LIST:
        tag_name = f"NiO_slab_{''.join(map(str, MILLER_INDICES))}_{n_layers}L"
        os.makedirs(tag_name, exist_ok=True)

        slab = surface(
            bulk,
            indices=MILLER_INDICES,
            layers=n_layers,
            vacuum=VACUUM_A,
            periodic=True,
        )
        slab.center(vacuum=VACUUM_A, axis=2)

        z = slab.positions[:, 2]
        z_cut = z.min() + FIX_BOTTOM_FRACTION * (z.max() - z.min())
        fixed_mask = z < z_cut
        slab.set_constraint(FixAtoms(mask=fixed_mask))

        n_ni1 = int(np.sum(slab.get_tags() == 1))
        n_ni2 = int(np.sum(slab.get_tags() == 2))
        n_o = int(np.sum(slab.get_tags() == 0))
        actual_vac = slab.cell[2, 2] - (z.max() - z.min())

        print(
            f"[{tag_name}] Ni1:Ni2:O = {n_ni1}:{n_ni2}:{n_o} | "
            f"Fixed = {int(np.sum(fixed_mask))}/{len(slab)} | Vac = {actual_vac:.2f} A"
        )

        write(os.path.join(tag_name, f"{tag_name}_init.cif"), slab)
        write_qe_slab_input(
            os.path.join(tag_name, f"{tag_name}.relax.in"),
            slab,
            fixed_mask,
            tag_name,
            SLAB_KPTS,
            PSEUDO_DIR,
        )
