cd ../result
awk -F, 'NR==FNR {a[$1","$2]; next} ($1","$2) in a' adapter_rows_sequence.csv adapter_rows_dbscan.csv > adapter.csv 
python ../process/convert.py


