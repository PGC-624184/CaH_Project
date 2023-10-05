#!/bin/bash

#factor
factor=45.56335252907954
start=422.673
# Create an output file
output_file="../Ca_H2_output.csv"
output_file_core="../Ca_H2_LineCore.csv"

datapath="../data/CaH2"
# Loop through the list of input files
for input_file in "$datapath"/*.csv; do
    # Use awk to extract the specified lines and calculate the result
    value_x=$(awk 'NR==7 {print $2}' "$input_file")
    value_y=$(awk 'NR==8 {print $2}' "$input_file")
    value_z=$(awk 'NR==9 {print $2}' "$input_file")
    result_x=$(echo $factor/$value_x - $start | bc -l )
    result_y=$(echo $factor/$value_y - $start | bc -l )
    result_z=$(echo $factor/$value_z - $start | bc -l )
    # Append the result to the output file
    filebase=$(basename "$input_file")
    number=$(echo "${filebase:25}")
    number_1=${number%.*}
    echo "$number_1, $result_x, $result_y, $result_z " >> "$output_file"
done