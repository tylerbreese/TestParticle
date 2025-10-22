#!/bin/bash
# Usage:
#   ./batch_run_matlab.sh [directory] [parallel_jobs]
#
# Examples:
#   ./batch_run_matlab.sh              # runs all .in files in current dir, sequentially
#   ./batch_run_matlab.sh inputs 4     # runs all .in files in 'inputs/' using 4 parallel jobs

set -e

# --- Arguments ---
INPUT_DIR="${1:-.}"               # Default: current directory
MAX_JOBS="${2:-1}"                # Default: run 1 job at a time (no parallelism)

# --- Collect .in files ---
IN_FILES=("$INPUT_DIR"/*.in)
NUM_FILES=${#IN_FILES[@]}

if [ "$NUM_FILES" -eq 0 ]; then
    echo "No .in files found in $INPUT_DIR"
    exit 1
fi

echo "Found $NUM_FILES input files in '$INPUT_DIR'"
echo "Running up to $MAX_JOBS MATLAB jobs in parallel..."
echo

# --- Function to run one MATLAB job ---
run_matlab_job() {
    local INPUT_FILE="$1"
    echo "[`date +'%H:%M:%S'`] Starting job for $INPUT_FILE"
    matlab -batch "input_file='${INPUT_FILE}'; run('run_helium.m')" > "${INPUT_FILE%.in}.log" 2>&1
    echo "[`date +'%H:%M:%S'`] Finished job for $INPUT_FILE"
}

# --- Run sequentially or in parallel ---
if [ "$MAX_JOBS" -le 1 ]; then
    # Sequential run
    for f in "${IN_FILES[@]}"; do
        run_matlab_job "$f"
    done
else
    # Parallel run
    job_count=0
    for f in "${IN_FILES[@]}"; do
        run_matlab_job "$f" &
        ((job_count++))

        # Limit number of concurrent jobs
        if (( job_count % MAX_JOBS == 0 )); then
            wait
        fi
    done
    wait
fi

echo
echo "✅ All MATLAB runs completed!"
