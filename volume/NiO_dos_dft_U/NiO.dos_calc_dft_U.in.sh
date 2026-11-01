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
DOS_BIN=dos.x
PROJ_BIN=projwfc.x

echo "PSEUDO DIR: $PSEUDO_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"

run_scf() {
    local INPUT_FILE="NiO.scf_dft_U.in"
    local OUTPUT_FILE="NiO.scf_dft_U.out"

    echo "Started scf calculation"

    # Using cat heredoc to generate the file content
    cat <<EOF >$INPUT_FILE
&CONTROL
    calculation = 'scf'
    prefix      = 'NiO'
    restart_mode = 'from_scratch'
    outdir      = '$OUTPUT_DIR'
    pseudo_dir  = '$PSEUDO_DIR'
    tprnfor     = .true.
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
/
&ELECTRONS
    mixing_beta = 0.1
    electron_maxstep = 300
    mixing_ndim = 14
    conv_thr    = 1d-9
/
&IONS
/
&CELL
    cell_dofree='ibrav'
/
ATOMIC_SPECIES
Ni1  58.6934  Ni.upf
Ni2  58.6934  Ni.upf
O   15.999   O.upf

ATOMIC_POSITIONS (crystal)
Ni1              0.0000000000        0.0000000000       -0.0000000000
Ni2             -0.5000000000        1.5000000000       -0.5000000000
O                0.7500000023       -0.2499999977       -0.2499999977
O               -0.7500000023        0.2499999977        0.2499999977

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

run_nscf() {
    local INPUT_FILE="NiO.nscf_dft_U.in"
    local OUTPUT_FILE="NiO.nscf_dft_U.out"

    echo "Started nscf calculation"

    # Using cat heredoc to generate the file content
    cat <<EOF >$INPUT_FILE
&CONTROL
    calculation = 'nscf'
    prefix      = 'NiO'
    restart_mode = 'from_scratch'
    outdir      = '$OUTPUT_DIR'
    pseudo_dir  = '$PSEUDO_DIR'
    tprnfor     = .true.
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
    occupations = 'tetrahedra'
/
&ELECTRONS
    mixing_beta = 0.1
    electron_maxstep = 300
    mixing_ndim = 14
    conv_thr    = 1d-9
/
&IONS
/
&CELL
    cell_dofree='ibrav'
/
ATOMIC_SPECIES
Ni1  58.6934  Ni.upf
Ni2  58.6934  Ni.upf
O   15.999   O.upf

ATOMIC_POSITIONS (crystal)
Ni1              0.0000000000        0.0000000000       -0.0000000000
Ni2             -0.5000000000        1.5000000000       -0.5000000000
O                0.7500000023       -0.2499999977       -0.2499999977
O               -0.7500000023        0.2499999977        0.2499999977

K_POINTS automatic
12 12 12 0 0 0

HUBBARD (ortho-atomic)
U Ni1-3d 5
U Ni2-3d 5
U O-2p 6
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 28 "$QE_BIN" -nk 4 < $INPUT_FILE > $OUTPUT_FILE
}

run_dos() {
    local INPUT_FILE="NiO.dos_dft_U.in"
    local OUTPUT_FILE="NiO.dos_dft_U.out"

    echo "Started dos calculation"

    # Using cat heredoc to generate the file content
    cat <<EOF >$INPUT_FILE
    &DOS
    prefix='NiO'
    outdir='$OUTPUT_DIR'
    fildos='nio.dos_dft_U.dat'
    emin=-10
    emax=20
/
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 28 "$DOS_BIN" -nk 4 < $INPUT_FILE > $OUTPUT_FILE
}

run_projwfc() {
    local INPUT_FILE="NiO.profwfc.in"
    local OUTPUT_FILE="NiO.profwfc.out"
    echo "Started projwfc calculation for local magnetic vectors"
    cat <<EOF >$INPUT_FILE
&PROJWFC
    prefix  = 'NiO'
    outdir  = '$OUTPUT_DIR'
    ngauss  = 0
    degauss = 0.01
    DeltaE  = 0.05
    Emin    = -10.0
    Emax    = 20.0
/
EOF
    rm -f $OUTPUT_FILE
    mpirun -np 28 "$PROJ_BIN" -nk 4 < $INPUT_FILE > $OUTPUT_FILE
}

run_plotter_dos() {
    local INPUT_FILE="plotter.gnu"

    echo "Plotting calculation results"

    # Using cat heredoc to generate the file content
    cat <<'EOF' >$INPUT_FILE
# Gnuplot script file for plotting data in file "si.dos.dat"
# This file is called plot_dos.gnu
set terminal pngcairo enhanced font 'Verdana,10'
set output 'NiO_dos_dft_U.png'
set grid
set title "Density of states (DOS) of NiO crystal"
set xlabel "Energy (eV)"
set ylabel "DOS"
set format y "%.4f"
plot    "nio.dos_dft_U.dat" using 1:2 title 'Total DOS' with line,\
	"nio.dos_dft_U.dat" using 1:3 title 'Integrated DOS' with line
EOF

    gnuplot ${INPUT_FILE}
}

run_plotter_projwfc() {
    local INPUT_FILE="plotter.gnu"
    echo "Plotting localized m_z magnetic moments"
    cat <<'EOF' >$INPUT_FILE
    set terminal pngcairo enhanced font 'Verdana,10' size 800,600
    set output 'NiO_pdos_mz.png'
    set grid
    set title "Projected Local Magnetization (m_z) of Ni 3d in NiO"
    set xlabel "Energy (eV)"
    set ylabel "PDOS (States / eV)"
    set format y "%.2f"
    
    plot "NiO.pdos_atm#1(Ni1)_wfc#4(d_j2.5)" using 1:5 title 'Ni1 3d j=5/2 (m_z)' with line lw 2 lc rgb 'red', \
         "NiO.pdos_atm#1(Ni1)_wfc#5(d_j1.5)" using 1:5 title 'Ni1 3d j=3/2 (m_z)' with line lw 2 lc rgb 'dark-red', \
         "NiO.pdos_atm#2(Ni2)_wfc#4(d_j2.5)" using 1:5 title 'Ni2 3d j=5/2 (m_z)' with line lw 2 lc rgb 'blue', \
         "NiO.pdos_atm#2(Ni2)_wfc#5(d_j1.5)" using 1:5 title 'Ni2 3d j=3/2 (m_z)' with line lw 2 lc rgb 'dark-blue'
EOF
    gnuplot ${INPUT_FILE}
}

run_scf
run_nscf
run_dos
run_projwfc
run_plotter_dos
run_plotter_projwfc
