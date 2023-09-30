#!/bin/bash

#factor
factor=45.56335252907954 
# Create an output file
output_file="../Ca_H2_output.csv"

datapath="../data/CaH2"
# Loop through the list of input files
for input_file in "$datapath"/*.csv; do
    # Use awk to extract the specified lines and calculate the result
    value_7=$(awk 'NR==7 {print $2}' "$input_file")
    value_9=$(awk 'NR==9 {print $2}' "$input_file")
    result=$(echo $factor/$value_7 - $factor/$value_9 | bc -l )
    # Append the result to the output file
    filebase=$(basename "$input_file")
    number=$(echo -n "${filebase:25:-4}")
    echo "$number, $result " >> "$output_file"
done