#!/bin/bash

set -xe

PSEUDO_DIR="$(pwd)/../../../pot"
OUTPUT_DIR="$(pwd)/out"

QE_DIR="$HOME/Desktop/src/q-e-qe-7.5/bin"

echo "PSEUDO DIR: $PSEUDO_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"

run_pw() {
    local magn=$1
    local u_val=$2
    local magn_1=$(echo "$magn * -1" | bc -l)
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
    starting_magnetization(1) = $magn
    starting_magnetization(2) = $magn_1 
    tot_magnetization = 0.0
    ecutwfc = 55
    ecutrho = 550
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
U Ni1-3d $u_val
U Ni2-3d $u_val
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 5 "$QE_DIR/pw.x" -inp $INPUT_FILE > $OUTPUT_FILE
}

for u in $(seq 3.0 0.5 8.0); do
    echo "Processing Hubbard U = $u"
    for m in $(seq 0.1 0.05 0.7); do
        echo "  Running Magn = $m..."
        run_pw "$m" "$u"
        ENERGY=$(grep "!" NiO.scf.out | tail -n 1 | awk '{print $5}')
        if [ -n "$ENERGY" ]; then
            echo "$m $u $ENERGY" >> NiO_heatmap.data
        fi
    done
    echo "" >> NiO_heatmap.data
done

# Gnuplot Heatmap Generation
echo "Generating Heatmap..."

cat <<EOF | gnuplot
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'energy_heatmap.png'

set title 'Total Energy Surface: Magnetization vs Hubbard U'
set xlabel 'Starting Magnetization'
set ylabel 'Hubbard U (eV)'

set view map
set pm3d at b
set palette rgbformulae 33,13,10  # Rainbow palette

# Adjust colorbar label
set cblabel "Total Energy (Ry)"

# Use splot for 3D data rendered as a map
splot 'NiO_heatmap.data' using 1:2:3 with pm3d title ''
EOF

echo "Heatmap saved to energy_heatmap.png"
