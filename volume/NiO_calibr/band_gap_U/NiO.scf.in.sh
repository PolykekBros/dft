#!/bin/bash

set -xe

PSEUDO_DIR="$(pwd)/../../../pot"
OUTPUT_DIR="$(pwd)/out"

QE_DIR="$HOME/Desktop/src/q-e-qe-7.5/bin"

echo "PSEUDO DIR: $PSEUDO_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"

run_pw() {
    local u_val=$1
    # Bash requires explicit arithmetic syntax for math
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
    starting_magnetization(1) =  0.2
    starting_magnetization(2) = -0.2
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
U Ni1-3d 4.6
U Ni2-3d 4.6
U O-2p $u_val
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 5 "$QE_DIR/pw.x" -inp $INPUT_FILE > $OUTPUT_FILE
}

> NiO_gap.data

for u in $(seq 1.0 0.25 6.0); do
    echo "Processing Hubbard U = $u"
    
    run_pw "$u"

    GAP_LINE=$(grep "highest occupied, lowest unoccupied level (ev):" NiO.scf.out)

    if [ -n "$GAP_LINE" ]; then
        VAL1=$(echo $GAP_LINE | awk '{print $7}')
        VAL2=$(echo $GAP_LINE | awk '{print $8}')
    
        GAP=$(echo "$VAL2 - $VAL1" | bc -l | sed 's/-/ /' | awk '{print $1}')
        
        echo "U = $u, Gap = $GAP eV"
        echo "$u $GAP" >> NiO_gap.data
    else
        echo "U = $u, No gap found (check if system is metallic)"
        echo "$u 0.0" >> NiO_gap.data
    fi
done

# Gnuplot Line Graph Generation
echo "Generating Band Gap Plot..."

cat <<EOF | gnuplot
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'NiO_band_gap.png'

set title 'NiO Band Gap vs. Hubbard U'
set xlabel 'Hubbard U (eV)'
set ylabel 'Band Gap (eV)'

set grid
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5

plot 'NiO_gap.data' with linespoints ls 1 title 'Fundamental Gap'
EOF
