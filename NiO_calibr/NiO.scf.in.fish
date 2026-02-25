#!/usr/bin/env fish

set -g fish_trace 1

set PSEUDO_DIR "$(pwd)/../pot"
set OUTPUT_DIR "$(pwd)/out"

echo "PSEUDO DIR: $PSEUDO_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"

function run_pw
    set ecutwfc $argv[1]
    set ecutrho (math "$ecutwfc * 10")
    set INPUT_FILE "NiO.scf.in"
    set OUTPUT_FILE "NiO.scf.out"

    echo "
&CONTROL
    calculation = 'scf'
    prefix      = 'NiO'
    restart_mode = 'from_scratch'
    outdir      = '$OUTPUT_DIR'
    pseudo_dir  = '$PSEUDO_DIR'
    tprnfor     = .true.
/
&SYSTEM
    ibrav = 0
    nat   = 12
    ntyp  = 2
    nspin = 1
    ecutwfc = $ecutwfc
    ecutrho = $ecutrho
    nbnd  = 80
/
&ELECTRONS
    mixing_beta = 0.3
    conv_thr    = 1d-10
    electron_maxstep = 1000
/
ATOMIC_SPECIES
Ni  58.6934  Ni.pbesol-n-rrkjus_psl.0.1.UPF
O   15.999   O.pbesol-n-rrkjus_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
Ni  0.00000000  0.00000000  0.00000000 
Ni  0.33333333  0.66666700  0.12880400 
Ni  0.66666700  0.33333300  0.25760900 
Ni  0.00000000  0.00000000  0.19320600 
Ni  0.33333333  0.66666700  0.32201100 
Ni  0.66666700  0.33333300  0.06440200 
O   0.00000000  0.00000000  0.28981000 
O   0.33333333  0.66666700  0.03220100 
O   0.66666700  0.33333300  0.16100500 
O   0.00000000  0.00000000  0.09660300 
O   0.33333333  0.66666700  0.22540800 
O   0.66666700  0.33333300  0.35421200 

K_POINTS automatic
4 4 1 0 0 0

CELL_PARAMETERS angstrom
 2.954859   0.000000   0.000000
-1.477429   2.558990   0.000000
 0.000000   0.000000  37.040310

HUBBARD (ortho-atomic)
U Ni-3d 4.3
" >$INPUT_FILE

    rm -f $OUTPUT_FILE
    mpirun -n (math (nproc) / 2) pw.x -inp $INPUT_FILE >$OUTPUT_FILE
    
    set energy (grep "!" $OUTPUT_FILE | grep "total energy" | tail -n 1 | awk '{print $4}')
    echo "ecutwfc: $ecutwfc, Total Energy: $energy Ry"
    
    # Save to data file
    echo "$ecutwfc $energy" >> NiO.scf.data
end

# Initialize/Clear the data file
echo "# ecutwfc TotalEnergy" > NiO.scf.data

# Loop through the variants
for val in 50 55 60 65 70
    run_pw $val
end

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
