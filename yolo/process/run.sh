#!/bin/bash
python run.py  --bam_filename "/root/autodl-tmp/adapter/H_h/H_m_no_qualities_new.bam" \
               --bamview_file "/root/autodl-tmp/adapter/H_h/bamview-new.csv"\
               --process_files "/root/autodl-tmp/pmf/pmf_dev2/data/detection_results.csv" \
               --similarity 0.7 \
               --num_threads 1 \
               --eps 100 \
               --min_samples 2 \
               --k_size 8\
               --images_dir "/root/autodl-tmp/adapter/H_h/images_m" \
               --output_images_dir "/root/autodl-tmp/pmf/pmf_dev2/output_images" \
               --output_detection_csv_path "/root/autodl-tmp/pmf/pmf_dev2/data/detection_results.csv" \
               --adapter_output_dir "/root/autodl-tmp/pmf/pmf_dev2/result"\
               --is_inference 1
