#!/bin/bash
# Process little.csv and divide all the identified targets according to contig.
cd ../data
rm -rf ptg*
cd ../process
python processImageName.py
# Process according to the processed files to obtain the adapters contained in each contig.
python all_in_one.py  --bam_filename "../data/SRR10238607_mdbg.bam" --similarity 0.7 --num_threads 2 --eps 100 --min_samples 2 --k_size 8
# merge adapter

bash merge_sequence.sh
bash merge_cluster.sh
bash merge_all.sh


    



