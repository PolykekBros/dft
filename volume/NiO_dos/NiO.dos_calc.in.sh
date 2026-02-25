#!/bin/bash

set -xe

PSEUDO_DIR="$(pwd)/../../pot"
OUTPUT_DIR="$(pwd)/out"

QE_DIR="$HOME/Desktop/src/q-e-qe-7.5/bin"

echo "PSEUDO DIR: $PSEUDO_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"

run_scf() {
    local ecutwfc=$1
    # Bash requires explicit arithmetic syntax for math
    local ecutrho=$(($ecutwfc * 10))
    local INPUT_FILE="NiO.scf.in"
    local OUTPUT_FILE="NiO.scf.out"

    echo "Started scf calculation"

    # Using cat heredoc to generate the file content
    cat <<EOF >$INPUT_FILE
&CONTROL
    calculation = 'scf'
    prefix      = 'NiO'
    restart_mode = 'from_scratch'
    outdir      = '$OUTPUT_DIR'
    pseudo_dir  = '$PSEUDO_DIR'
/
&SYSTEM
    ibrav = 2
    celldm(1) = 7.837
    nat   = 2
    ntyp  = 2
    ecutwfc = $ecutwfc
    ecutrho = $ecutrho
    nbnd = 36,
/
&ELECTRONS
    mixing_beta = 0.3
    conv_thr    = 1d-9
/
ATOMIC_SPECIES
Ni  58.6934  Ni.pbesol-n-rrkjus_psl.0.1.UPF
O   15.999   O.pbesol-n-rrkjus_psl.1.0.0.UPF

ATOMIC_POSITIONS (alat)
Ni  0.00000000  0.00000000  0.00000000 0 0 0
O   0.50000000  0.50000000  0.50000000 0 0 0

K_POINTS automatic
4 4 4 0 0 0
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 5 "$QE_DIR/pw.x" -inp $INPUT_FILE > $OUTPUT_FILE
}

run_nscf() {
    local ecutwfc=$1
    # Bash requires explicit arithmetic syntax for math
    local ecutrho=$(($ecutwfc * 10))
    local INPUT_FILE="NiO.nscf.in"
    local OUTPUT_FILE="NiO.nscf.out"

    echo "Started nscf calculation"

    # Using cat heredoc to generate the file content
    cat <<EOF >$INPUT_FILE
&CONTROL
    calculation = 'nscf'
    prefix      = 'NiO'
    restart_mode = 'from_scratch'
    outdir      = '$OUTPUT_DIR'
    pseudo_dir  = '$PSEUDO_DIR'
/
&SYSTEM
    ibrav = 2
    celldm(1) = 7.837
    nat   = 2
    ntyp  = 2
    ecutwfc = $ecutwfc
    ecutrho = $ecutrho
    nbnd = 36,
    occupations = 'tetrahedra'
/
&ELECTRONS
    mixing_beta = 0.3
    conv_thr    = 1d-9
/
ATOMIC_SPECIES
Ni  58.6934  Ni.pbesol-n-rrkjus_psl.0.1.UPF
O   15.999   O.pbesol-n-rrkjus_psl.1.0.0.UPF

ATOMIC_POSITIONS (alat)
Ni  0.00000000  0.00000000  0.00000000 0 0 0
O   0.50000000  0.50000000  0.50000000 0 0 0

K_POINTS automatic
12 12 12 0 0 0
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 5 "$QE_DIR/pw.x" -inp $INPUT_FILE > $OUTPUT_FILE
}

run_dos() {
    local INPUT_FILE="NiO.dos.in"
    local OUTPUT_FILE="NiO.dos.out"

    echo "Started dos calculation"

    # Using cat heredoc to generate the file content
    cat <<EOF >$INPUT_FILE
    &DOS
    prefix='NiO'
    outdir='$OUTPUT_DIR'
    fildos='nio.dos.dat'
    emin=-10
    emax=20
/
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 5 "$QE_DIR/dos.x" -inp $INPUT_FILE > $OUTPUT_FILE
}

run_plotter() {
    local INPUT_FILE="plotter.gnu"

    echo "Plotting calculation results"

    # Using cat heredoc to generate the file content
    cat <<EOF >$INPUT_FILE
# Gnuplot script file for plotting data in file "si.dos.dat"
# This file is called plot_dos.gnu
set terminal pngcairo enhanced font 'Verdana,10'
set output 'NiO_dos.png'
set grid
set title "Density of states (DOS) of NiO crystal"
set xlabel "Energy (eV)"
set ylabel "DOS"
set format y "%.4f"
plot    "nio.dos.dat" using 1:2 with linespoints title 'DOS'
EOF

    gnuplot ${INPUT_FILE}
}

run_scf 55
run_nscf 55
run_dos
run_plotter
