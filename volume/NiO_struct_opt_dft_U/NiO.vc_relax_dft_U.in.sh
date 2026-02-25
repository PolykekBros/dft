#!/bin/bash

set -xe

PSEUDO_DIR="$(pwd)/../../pot"
OUTPUT_DIR="$(pwd)/out"

QE_DIR="$HOME/Desktop/src/q-e-qe-7.5/bin"

echo "PSEUDO DIR: $PSEUDO_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"

run_vc_relax() {
    local ecutwfc=$1
    # Bash requires explicit arithmetic syntax for math
    local ecutrho=$(($ecutwfc * 10))
    local INPUT_FILE="NiO.vc_relax_dft_U.in"
    local OUTPUT_FILE="NiO.vc_relax_dft_U.out"

    # Using cat heredoc to generate the file content
    cat <<EOF >$INPUT_FILE
&CONTROL
    calculation = 'vc-relax'
    prefix      = 'NiO'
    restart_mode = 'from_scratch'
    outdir      = '$OUTPUT_DIR'
    pseudo_dir  = '$PSEUDO_DIR'
    tprnfor     = .true.
    etot_conv_thr = 1e-5
    forc_conv_thr = 1e-4 
/
&SYSTEM
    ibrav = 2
    celldm(1) = 7.9
    nat   = 4
    ntyp  = 3
    nspin = 2
    starting_magnetization(1) =  0.5
    starting_magnetization(2) = -0.5
    tot_magnetization = 0.0
    ecutwfc = $ecutwfc
    ecutrho = $ecutrho
/
&ELECTRONS
    mixing_beta = 0.3
    conv_thr    = 1d-9
    electron_maxstep = 1000
/
&IONS
/
&CELL
    cell_dofree='all'
/
ATOMIC_SPECIES
 Ni1  58.6934  Ni.pbesol-n-rrkjus_psl.0.1.UPF
 Ni2  58.6934  Ni.pbesol-n-rrkjus_psl.0.1.UPF
 O    15.999   O.pbesol-n-rrkjus_psl.1.0.0.UPF

ATOMIC_POSITIONS (alat)
 Ni1  0.0  0.0  0.0 0 0 0
 Ni2  0.5  0.5  0.5 0 0 0
 O    0.25 0.25 0.25 0 0 0
 O    0.75 0.75 0.75 0 0 0 

K_POINTS (automatic)
4 4 4 0 0 0

HUBBARD (ortho-atomic)
U Ni1-3d 4.6
U Ni2-3d 4.6
EOF

    rm -f $OUTPUT_FILE
    # Using nproc to calculate processes, similar to the fish script
    mpirun -np 5 "$QE_DIR/pw.x" -inp $INPUT_FILE > $OUTPUT_FILE
}

run_vc_relax 55
