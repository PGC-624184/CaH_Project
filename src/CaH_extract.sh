#!/bin/bash

#factor
factor=45.56335252907954
start=422.673
# Create an output file
output_file="../Ca_H_output.csv"

datapath="../data/CaH"
# Loop through the list of input files
for input_file in "$datapath"/*.csv; do
    # Use awk to extract the specified lines and calculate the result
    value_x=$(awk 'NR==15 {print $2}' "$input_file")
    value_y=$(awk 'NR==16 {print $2}' "$input_file")
    value_z=$(awk 'NR==17 {print $2}' "$input_file")
    result_x=$(echo $factor/$value_x - $start | bc -l )
    result_y=$(echo $factor/$value_y - $start | bc -l )
    result_z=$(echo $factor/$value_z - $start | bc -l )
    filebase=$(basename "$input_file")
    number=$(echo -n "${filebase:24}")
    number_1=${number%.*}
    # Append the result to the output file
    echo "$number_1, $result_x, $result_y, $result_z " >> "$output_file"
done