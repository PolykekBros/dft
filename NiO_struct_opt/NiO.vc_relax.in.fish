#!/usr/bin/env fish

set -g fish_trace 1

set PSEUDO_DIR "$(pwd)/../pot"
set OUTPUT_DIR "$(pwd)/out"

echo "PSEUDO DIR: $PSEUDO_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"

function run_vc_relax
    set ecutwfc $argv[1]
    set ecutrho (math "$ecutwfc * 10")
    set INPUT_FILE "NiO.vc_relax.in"
    set OUTPUT_FILE "NiO.vc_relax.out"

    echo "
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
    conv_thr    = 1d-6
    electron_maxstep = 1000
/
&IONS
/
&CELL
    cell_dofree='ibrav'
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
4 4 4 0 0 0

CELL_PARAMETERS angstrom
        0.8616097902        -0.4974526223         2.7823301882
       -1.5409629259         2.4342380536        -0.6566327261
       -12.201827083        -7.0447001392         2.5190429035
" >$INPUT_FILE

    rm -f $OUTPUT_FILE
    mpirun -np (math (nproc) / 2) pw.x -inp $INPUT_FILE >$OUTPUT_FILE
end

run_vc_relax 55
