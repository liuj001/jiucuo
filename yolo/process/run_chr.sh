#!/bin/bash
python run.py  --bam_filename "/root/autodl-tmp/HG002/adapter/chr21_adapter.bam" \
               --bamview_file "/root/autodl-tmp/HG002/adapter/bamview-new.csv"\
               --images_dir "/root/autodl-tmp/HG002/adapter/images_0106" \
               --similarity 0.7 \
               --num_threads 1 \
               --eps 100 \
               --min_samples 2\
               --k_size 8\
               --adapter_output_dir "../result_chr_0106_images"\
               --is_inference 1

