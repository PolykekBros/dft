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
    nbnd = 32,
/
&ELECTRONS
    mixing_beta = 0.3
    conv_thr    = 1d-9
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
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 5 "$QE_DIR/pw.x" -inp $INPUT_FILE > $OUTPUT_FILE
}

run_nscf() {
    local ecutwfc=$1
    # Bash requires explicit arithmetic syntax for math
    local ecutrho=$(($ecutwfc * 10))
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
    nbnd = 32,
    occupations = 'tetrahedra'
/
&ELECTRONS
    mixing_beta = 0.3
    conv_thr    = 1d-9
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
12 12 12 0 0 0
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 5 "$QE_DIR/pw.x" -inp $INPUT_FILE > $OUTPUT_FILE
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
    mpirun -np 5 "$QE_DIR/dos.x" -inp $INPUT_FILE > $OUTPUT_FILE
}

run_plotter() {
    local INPUT_FILE="plotter.gnu"

    echo "Plotting calculation results"

    # Using cat heredoc to generate the file content
    cat <<EOF >$INPUT_FILE
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Density of states (DOS) of GaAs with DFT+U"
set xlabel "Energy (eV)"
set ylabel "DOS"
set arrow 1 from 15.691,-2.5 to 15.691,2.5 nohead ls 10 dt 2
#set xr [5:25]
#set yr [0:325]
plot    "nio.dos_dft_U.dat" using 1:2 title 'DOS of spin-up' with line,\
	"nio.dos_dft_U.dat" using 1:(-$3) title 'DOS of spin-down' with line
pause -1 "Hit any key to continue\n"    #so that the code doesn't exit automatically
EOF

    gnuplot ${INPUT_FILE}
}

run_scf 55
run_nscf 55
run_dos
run_plotter
