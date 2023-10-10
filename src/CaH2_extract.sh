#!/bin/bash

# Create an output file
output_file="../data/Ca_H2_test.csv"

datapath="../data/CaH2"
# Loop through the list of input files
for input_file in "$datapath"/*.csv; do
    # Use awk to extract the specified lines and calculate the result
    # Oscillator strength over 0.05, wavelength between 300-500nm
    values=$(awk '$3>=0.05 && $2>=0.09 && $2<=0.15 {print $2}' "$input_file" | paste -s -)
    # Append the result to the output file
    filebase=$(basename "$input_file")
    number=$(echo "${filebase:25}")
    number_1=${number%.*}
    echo  -e "$number_1 \t $values" >> "$output_file"
done