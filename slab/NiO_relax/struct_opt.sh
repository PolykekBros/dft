#!/bin/bash

#SBATCH --job-name=NiO_slab_scan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --output=NiO_slab.out

# --- Environment Setup ---
module purge
module load mpi/impi/2021.4.0
module load espresso/7.4.1
module load post/gnuplot
set -xe

unset I_MPI_PMI_LIBRARY
export I_MPI_FABRICS=shm:ofi

QE_BIN=pw.x
SLAB_DIR="NiO_slab_100_4L"
INPUT_FILE="${SLAB_DIR}.relax.in"
OUTPUT_FILE="${SLAB_DIR}.relax.out"

# --- Generate Slab Inputs ---
uv run slab_gen.py

# --- Run 4L Slab Relaxation ---
cd "$SLAB_DIR"
mkdir -p out
rm -f "$OUTPUT_FILE"

mpirun -np 28 "$QE_BIN" -nk 4 < "$INPUT_FILE" > "$OUTPUT_FILE"
cd ..

# --- Post-Relaxation Verification ---
grep -q "JOB DONE." "${SLAB_DIR}/${OUTPUT_FILE}"

uv run python -c "
from slab_gen import verify_slab_convergence, BULK_OUT
verify_slab_convergence(BULK_OUT, ['${SLAB_DIR}/${OUTPUT_FILE}'])
"
