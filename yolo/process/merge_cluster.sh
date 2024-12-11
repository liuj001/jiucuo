#!/bin/bash
#  This script merges all output results from contigs in the ../data directory and outputs them as adapter_rows_dbscan.csv
cd ../data

folder_path="."


output_file="$(pwd)/adapter_rows_dbscan.csv"
touch "$output_file"


first_file=$(find "$folder_path" -maxdepth 1 -type f -name '*second_part.csv' | sort | head -n 1)


if [ -n "$first_file" ]; then
    echo "--------------------"
    head -n 1 "$first_file" > "$output_file"
   
fi


for file in "$folder_path"/*second_part.csv; do
    if [ -e "$file" ] ; then
        
        tail -n +2 "$file" >> "$output_file"
        
    fi
done

echo "saved file : $output_file"
mv $output_file ../result



