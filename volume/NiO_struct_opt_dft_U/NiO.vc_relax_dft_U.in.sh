#!/bin/bash

#SBATCH --job-name=NiO_U_scan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --output=NiO_scan.out

## --- Environment Setup ---
module purge
module load mpi/impi/2021.4.0
module load espresso/7.4.1
module load post/gnuplot
set -xe

unset I_MPI_PMI_LIBRARY
export I_MPI_FABRICS=shm:ofi

PSEUDO_DIR="$(pwd)/../../pot"
OUTPUT_DIR="$(pwd)/out"
QE_BIN=pw.x

mkdir -p "$OUTPUT_DIR"

echo "PSEUDO DIR: $PSEUDO_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"

run_vc_relax() {
    local INPUT_FILE="NiO.vc_relax_dft_U.in"
    local OUTPUT_FILE="NiO.vc_relax_dft_U.out"

    # Using cat heredoc to generate the file content
    cat <<EOF >$INPUT_FILE
&CONTROL
    calculation = 'vc-relax'
    prefix      = 'NiO'
    restart_mode = 'from_scratch'
    outdir      = '$OUTPUT_DIR'
    pseudo_dir  = '$PSEUDO_DIR'
    tprnfor     = .true.
    etot_conv_thr = 1e-5
    forc_conv_thr = 1e-4 
/
&SYSTEM
    ibrav = 5
    celldm(1) = 9.711
    celldm(4) = 0.8333
    nat   = 4
    ntyp  = 3
    noncolin = .true.
    lspinorb = .true.
    starting_magnetization(1) = 0.5
    starting_magnetization(2) = 0.5
    angle1(1) = 0.0
    angle1(2) = 180.0
    angle2(1) = 0.0
    angle2(2) = 0.0

    ecutwfc = 65
    ecutrho = 260
    nbnd = 80

    occupations = 'smearing'
    smearing = 'cold'
    degauss = 0.01
/
&ELECTRONS
    mixing_beta = 0.1
    electron_maxstep = 300
    mixing_ndim = 14
    conv_thr    = 1d-9
/
&IONS
    ion_dynamics='bfgs'
/
&CELL
    cell_dofree='all'
    cell_dynamics='bfgs'
/
ATOMIC_SPECIES
 Ni1  58.6934  Ni.upf
 Ni2  58.6934  Ni.upf
 O    15.999   O.upf

ATOMIC_POSITIONS (crystal)
    Ni1  0.0000000000   0.0000000000   0.0000000000
    Ni2 -0.5000000000   1.5000000000  -0.5000000000
    O    0.7500000000  -0.2500000000  -0.2500000000
    O   -0.7500000000   0.2500000000   0.2500000000

K_POINTS automatic
6 6 6 0 0 0

HUBBARD (ortho-atomic)
U Ni1-3d 5
U Ni2-3d 5
U O-2p 6
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 28 "$QE_BIN" -nk 4 < $INPUT_FILE > $OUTPUT_FILE
}

run_vc_relax

> NiO_structural.data

DENSITY=$(grep "density =" NiO.vc_relax_dft_U.out | tail -n 1 | awk '{print $3}')
ALAT_VAL=$(grep "CELL_PARAMETERS (alat=" NiO.vc_relax_dft_U.out | tail -n 1 | tr -d ')' | awk -F'=' '{print $2}')
RAW_PARAM=$(grep -A 1 "CELL_PARAMETERS (alat=" NiO.vc_relax_dft_U.out | tail -n 1 | awk '{print $1}')
CALC_PARAM_B=$(echo "scale=6; $RAW_PARAM / -0.5 * $ALAT_VAL" | bc -l)
CALC_PARAM_A=$(echo "scale=6; $CALC_PARAM_B * 0.529177" | bc -l)

echo "Density=$DENSITY | Calc_Param=$CALC_PARAM_A"
echo -e  "Density: $DENSITY\nLattice param (Bohr): $CALC_PARAM_B\nLattice param (Angstrom): $CALC_PARAM_A" >> NiO_structural.data

