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
BANDS_BIN=bands.x

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

run_bands() {
    local INPUT_FILE="NiO.bands_dft_U.in"
    local OUTPUT_FILE="NiO.bands_dft_U.out"

    echo "Started nscf calculation"
    cat <<EOF >$INPUT_FILE
&CONTROL
    calculation = 'bands'
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
ATOMIC_SPECIES
 Ni1  58.6934  Ni.upf
 Ni2  58.6934  Ni.upf
 O    15.999   O.upf

ATOMIC_POSITIONS (crystal)
Ni1              0.0000000000        0.0000000000       -0.0000000000
Ni2             -0.5000000000        1.5000000000       -0.5000000000
O                0.7500000023       -0.2499999977       -0.2499999977
O               -0.7500000023        0.2499999977        0.2499999977

! High-symmetry path using crystal coordinates (crystal_b)
! The last integer dictates the number of k-points generated along that segment
K_POINTS crystal_b
6
  0.500  0.000  0.500  20 ! F
  0.000  0.000  0.000  20 ! Gamma
  0.500  0.500  0.500  20 ! T
  0.375  0.375  0.750  20 ! K
  0.000  0.000  0.500  20 ! L
  0.000  0.000  0.000  0  ! Gamma

HUBBARD (ortho-atomic)
U Ni1-3d 5
U Ni2-3d 5
U O-2p 6
EOF
    
    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 28 "$QE_BIN" -nk 4 < $INPUT_FILE > $OUTPUT_FILE
}

run_bands_pp() {
    local INPUT_FILE="NiO.bands_dft_U.pp.in"
    local OUTPUT_FILE="NiO.bands_dft_U.pp.out"

    echo "Started nscf calculation"
    cat <<EOF >$INPUT_FILE
&BANDS
prefix  = 'NiO'
outdir  = '$OUTPUT_DIR'
filband = 'NiO_bands.dat'
lsym    = .true.
/
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 28 "$BANDS_BIN" -nk 4 < $INPUT_FILE > $OUTPUT_FILE
}

run_plotter() {
    local INPUT_FILE="plotter.gnu"
    FERMI_ENERGY=$(grep "the Fermi energy is" NiO.nscf_dft_U.out | awk '{print $5}')
    echo "Extracted Fermi Energy: $FERMI_ENERGY eV"
    X_COORDS=($(grep "high-symmetry point" NiO.bands_dft_U.pp.out | awk '{print $NF}'))
    echo "K-point X-coordinates: ${X_COORDS[@]}"

    echo "Plotting band structure"
    cat <<EOF >$INPUT_FILE
    set terminal pngcairo enhanced font 'Verdana,12' size 800,600
    set output 'NiO_bands_relativistic.png'
    
    set title "Relativistic Band Structure of NiO (DFT+U+SOC)"
    set ylabel "Energy - E_F (eV)"
    set yrange [-8:8]  # Focus specifically on the Mott gap region
    set xzeroaxis ls -1
    set grid y
    
    set xtics ("L" ${X_COORDS[0]}, "{/Symbol G}" ${X_COORDS[1]}, "X" ${X_COORDS[2]}, "W" ${X_COORDS[3]}, "{/Symbol G}" ${X_COORDS[4]})
    
    # Add vertical dashed lines at the high-symmetry points for clarity
    set arrow from ${X_COORDS[1]}, graph 0 to ${X_COORDS[1]}, graph 1 nohead lt 0
    set arrow from ${X_COORDS[2]}, graph 0 to ${X_COORDS[2]}, graph 1 nohead lt 0
    set arrow from ${X_COORDS[3]}, graph 0 to ${X_COORDS[3]}, graph 1 nohead lt 0
    
    # Plotting the pre-formatted GNU output from bands.x
    plot 'NiO_bands.dat.gnu' using 1:(\$2 - ${FERMI_ENERGY}) with lines lw 1.5 lc rgb 'dark-blue' notitle
EOF
    gnuplot ${INPUT_FILE}
}

run_scf
run_nscf
run_bands
run_bands_pp
run_plotter
