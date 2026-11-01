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

PSEUDO_DIR="$(pwd)/../../../pot"
OUTPUT_DIR="$(pwd)/out"
QE_BIN=pw.x

mkdir -p "$OUTPUT_DIR"

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
    ecutrho = 650
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
    Ni1  0.0000000000   0.0000000000   0.0000000000
    Ni2 -0.5000000000   1.5000000000  -0.5000000000
    O    0.7500000000  -0.2500000000  -0.2500000000
    O   -0.7500000000   0.2500000000   0.2500000000

K_POINTS automatic
6 6 6 0 0 0

HUBBARD (ortho-atomic)
U Ni1-3d 5
U Ni2-3d 5
U O-2p $u_val
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 28 "$QE_BIN" -nk 4 < $INPUT_FILE > $OUTPUT_FILE
}

> NiO_gap.data

for u in $(seq 3.0 0.5 8.0); do
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
