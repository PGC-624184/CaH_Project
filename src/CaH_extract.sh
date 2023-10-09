#!/bin/bash

# Create an output file
output_file="../Ca_H_test.csv"

datapath="../data/CaH"
# Loop through the list of input files
for input_file in "$datapath"/*.csv; do
    # Use awk to extract the specified lines and calculate the result
    # choose oscillator strength >=0.07
    # choose wavelengths between 506 nm aand 350nm
    values=$(awk '$3>=0.05 && $2>=0.09 && $2<=0.15 {print $2}' "$input_file" | paste -s -)
    filebase=$(basename "$input_file")
    number=$(echo -n "${filebase:24}")
    number_1=${number%.*}
    # Append the result to the output file
    echo -e "$number_1 \t $values " >> "$output_file"
done