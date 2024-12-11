#!/bin/bash
#  This script merges all output results from contigs in the ../result directory and outputs them as adapter_rows_sequence.csv
cd ../result
folder_path="."


output_file="$(pwd)/adapter_rows_sequence.csv"
touch "$output_file"

# find head.csv 
first_file=$(find "$folder_path" -maxdepth 1 -type f -name 'head.csv' | sort | head -n 1)


if [ -n "$first_file" ]; then
    echo "--------------------"
    head -n 1 "$first_file" > "$output_file"
    
fi

# Loop through CSV files starting with "adapter_rows" (excluding the first file) and append them to the merged file (skipping the header)
for file in "$folder_path"/*all_detections_to_pre_adapter.csv; do
    if [ -e "$file" ] ; then   
        tail -n +2 "$file" >> "$output_file"
    fi
done

echo "saved file : $output_file"
rm -rf *all_detections_to_pre_adapter.csv



