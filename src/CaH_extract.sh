#!/bin/bash

#factor
factor=45.56335252907954
start=4.0   
# Create an output file
output_file="../Ca_H_output.csv"

datapath="../data/CaH"
# Loop through the list of input files
for input_file in "$datapath"/*.csv; do
    # Use awk to extract the specified lines and calculate the result
    value_7=$(awk 'NR==15 {print $2}' "$input_file")
    value_9=$(awk 'NR==17 {print $2}' "$input_file")
    result=$(echo $factor/$value_7 - $factor/$value_9 | bc -l )
    filebase=$(basename "$input_file")
    number=$(echo -n "${filebase:24:-4}")
    # Append the result to the output file
    echo "$number, $result " >> "$output_file"
done