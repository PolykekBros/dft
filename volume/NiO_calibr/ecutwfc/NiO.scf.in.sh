#!/bin/bash

set -xe

PSEUDO_DIR="$(pwd)/../../../pot"
OUTPUT_DIR="$(pwd)/out"

QE_DIR="$HOME/Desktop/src/q-e-qe-7.5/bin"

echo "PSEUDO DIR: $PSEUDO_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"

run_pw() {
    local ecutwfc=$1
    # Bash requires explicit arithmetic syntax for math
    local ecutrho=$(($ecutwfc * 10))
    local INPUT_FILE="NiO.scf.in"
    local OUTPUT_FILE="NiO.scf.out"

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
    ibrav = 2
    celldm(1) = 7.837
    nat   = 4
    ntyp  = 3
    nspin = 2
    starting_magnetization(1) =  0.5
    starting_magnetization(2) = -0.5
    tot_magnetization = 0.0
    ecutwfc = $ecutwfc
    ecutrho = $ecutrho
    nbnd = 32
/
&ELECTRONS
    mixing_beta = 0.3
    conv_thr    = 1d-9
/
&IONS
/
&CELL
    cell_dofree='ibrav'
/
ATOMIC_SPECIES
Ni1  58.6934  Ni.pbesol-n-rrkjus_psl.0.1.UPF
Ni2  58.6934  Ni.pbesol-n-rrkjus_psl.0.1.UPF
O   15.999   O.pbesol-n-rrkjus_psl.1.0.0.UPF

ATOMIC_POSITIONS (alat)
 Ni1  0.0  0.0  0.0 0 0 0
 Ni2  0.5  0.5  0.5 0 0 0
 O    0.25 0.25 0.25 0 0 0
 O    0.75 0.75 0.75 0 0 0 

K_POINTS automatic
4 4 4 0 0 0

HUBBARD (ortho-atomic)
U Ni1-3d 4.6
U Ni2-3d 4.6
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 5 "$QE_DIR/pw.x" -inp $INPUT_FILE > $OUTPUT_FILE
}

> NiO.scf.data 
for val in $(seq 30 5 70); do
    echo "Running calculation for ecutwfc = $val..."
    run_pw "$val"
    ENERGY=$(grep "!" NiO.scf.out | tail -n 1 | awk '{print $5}')
    if [ -z "$ENERGY" ]; then
        echo "Error: Could not find energy for $val Ry. Check NiO.vc_relax.out"
    else
        echo "Final Energy for $val Ry: $ENERGY Ry"
        echo "$val $ENERGY" >> NiO.scf.data
    fi
done

# Gnuplot to visualize results
echo "Generating convergence plot..."

cat <<EOF | gnuplot
set terminal pngcairo enhanced font 'Verdana,10'
set output 'convergence.png'
set title 'Total Energy Convergence vs. ecutwfc'
set xlabel 'ecutwfc (Ry)'
set ylabel 'Total Energy (Ry)'
set grid
set format y "%.4f"
plot 'NiO.scf.data' with linespoints title 'Converged Energy'
EOF

echo "Plot saved to convergence.png"

