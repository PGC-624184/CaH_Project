#!/bin/bash

#factor
factor=45.56335252907954
start=4.0   
# Create an output file
output_file="../output.csv"

datapath="../data"
# Loop through the list of input files
for input_file in "$datapath"/*.csv; do
    # Use awk to extract the specified lines and calculate the result
    value_7=$(awk 'NR==7 {print $2}' "$input_file")
    value_9=$(awk 'NR==9 {print $2}' "$input_file")
    result=$(echo $factor/$value_7 - $factor/$value_9 | bc -l )
    # Append the result to the output file
    echo "$start, $result " >> "$output_file"
    start=$(echo $start + 0.25 | bc -l)
done