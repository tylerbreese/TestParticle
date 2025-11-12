#!/bin/bash
# Usage: ./run_matlab.sh inputfile.in

# Exit if any command fails
set -e

# Check for argument
if [ $# -lt 1 ]; then
    echo "Usage: $0 input_file"
    exit 1
fi

# Get the input filename from first argument
INPUT_FILE=$1

# Optional: check that the file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File '$INPUT_FILE' not found!"
    exit 1
fi

# Run MATLAB in batch mode
echo "Running MATLAB with input file: $INPUT_FILE"
matlab -batch "input_file='${INPUT_FILE}'; run('multi.m')"


# make exec
# chmod +x batch_run_matlab.sh

# Example run once
# ./run_matlab.sh input.in
# Example run several in parallel
# ./run_matlab.sh input1.in &
# ./run_matlab.sh input2.in &
# ./run_matlab.sh input3.in &
# wait
# To stop jobs pkill -f matlab
