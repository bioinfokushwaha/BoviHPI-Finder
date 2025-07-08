#!/bin/bash
set -e

echo "Running bovihpi pipeline..."

# Use tools from Exe folder
chmod +x Exe/*

# Example run of a provided tool (replace with your actual flow)
./Exe/aa_calc_0 input.fasta > output_aa0.txt
./Exe/cc_calc input.fasta > output_cc.txt

# Add full toolchain workflow here...
echo "Pipeline completed successfully."
