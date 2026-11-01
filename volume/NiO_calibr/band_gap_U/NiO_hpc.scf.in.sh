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
QE_BIN="pw.x" 

mkdir -p "$OUTPUT_DIR"

## --- The Runner Function ---
run_pw() {
    local u_val=$1
    local INPUT_FILE="NiO.scf.in"
    local OUTPUT_FILE="NiO.scf.out"

    # Modern QE 7.x Hubbard syntax uses the HUBBARD card
    cat <<EOF > "$INPUT_FILE"
&CONTROL
    calculation  = 'scf'
    prefix       = 'NiO'
    restart_mode = 'from_scratch'
    outdir       = '$OUTPUT_DIR'
    pseudo_dir   = '$PSEUDO_DIR'
    tprnfor      = .true.
/
&SYSTEM
    ibrav = 1
    celldm(1) = 16
    nspin = 2
    nat  = 16
    ntyp = 3
    starting_magnetization(1) =  0.2
    starting_magnetization(2) = -0.2
    tot_magnetization = 0.0
    ecutwfc = 55.0
    ecutrho = 550.0
    nbnd = 120
/
&ELECTRONS
    mixing_beta = 0.3
    conv_thr    = 1d-9
/
&IONS
/

ATOMIC_SPECIES
 Ni1  58.6934  Ni.pbesol-n-rrkjus_psl.0.1.UPF
 Ni2  58.6934  Ni.pbesol-n-rrkjus_psl.0.1.UPF
 O    15.999   O.pbesol-n-rrkjus_psl.1.0.0.UPF

ATOMIC_POSITIONS (alat)
Ni2  0.00000  0.00000  0.00000
Ni1  0.00000  0.00000  0.50000
Ni1  0.00000  0.50000  0.00000
Ni2  0.00000  0.50000  0.50000
Ni1  0.50000  0.00000  0.00000
Ni2  0.50000  0.00000  0.50000
Ni2  0.50000  0.50000  0.00000
Ni1  0.50000  0.50000  0.50000
O   0.25000  0.25000  0.25000
O   0.25000  0.25000  0.75000
O   0.25000  0.75000  0.25000
O   0.25000  0.75000  0.75000
O   0.75000  0.25000  0.25000
O   0.75000  0.25000  0.75000
O   0.75000  0.75000  0.25000
O   0.75000  0.75000  0.75000

K_POINTS automatic
 2 2 2 0 0 0

HUBBARD (ortho-atomic)
 U Ni1-3d $u_val
 U Ni2-3d $u_val
EOF

    # Clean previous output to avoid grep confusion
    rm -f "$OUTPUT_FILE"
    
    # Execute with srun to use all allocated HPC cores
    mpirun -np 28 "$QE_BIN" -nk 4 < "$INPUT_FILE" > "$OUTPUT_FILE"
}

## --- Main Loop ---
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

## --- Visualization ---
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
