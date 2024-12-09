cd ../result
awk -F, 'NR==FNR {a[$1","$2]; next} ($1","$2) in a' adapter_rows.csv adapter_rows_second.csv > adapter.csv
