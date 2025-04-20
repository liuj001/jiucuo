# -*- coding: utf-8 -*-
import os
import argparse
import pysam
from time import *
#from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor, as_completed
from threading import Thread
import subprocess
import sys
import logging
import torch
import multiprocessing
from tqdm import tqdm
from adapter_speed_final import adapter_pic
from bamview_new import bamview
from add_pos import add_pos

parser=argparse.ArgumentParser(description='runJiuCuo.py integrates the processes of error correction.')
parser.add_argument('-output',type=str,help='Output directory',required=True)
parser.add_argument('-threads',type=int,help='Number of threads during correction',default=8)
parser.add_argument('-error_correction',type=int,help='Error correction',default=1)
args=parser.parse_args()
outdir = args.output
thread = args.threads
errc = args.error_correction

bam = outdir+'/raw.bam'
bam_dir = outdir+'/bam'
txt_dir = outdir+'/txt'
# txt_add_dir = outdir+'/txt_add'
bcf_txt_dir = outdir+'/bcf_txt'
snp_c_dir = outdir+'/snp_candidate'
adapter_dir = outdir+'/adapter'
adapter_out_dir = outdir+'/adapter_out'
bamview_txt_dir = outdir+'/bamview_txt'
bamview_csv_dir = outdir+'/bamview_csv'
snp_dir = outdir+'/snp'
ec_dir = outdir+'/ec_txt'
c_reads_dir = outdir+'/c_reads'
log_file = outdir+'/TOOLS_LOG.log'

def run_adapter(chr):
    chr = chr.rstrip('.bam')
    bamview_txt_cmd = "samtools view %s/%s.bam > %s/%s_bamview.txt"%(bam_dir, chr, bamview_txt_dir, chr)
    with open(log_file, "a") as log:
        bamview_txt_p=subprocess.Popen(bamview_txt_cmd,shell=True, stdout=log, stderr=log)
        bamview_txt_code=bamview_txt_p.wait()  #等待子进程结束，并返回状态码;
    csv_f=bamview_csv_dir+'/'+chr+'.csv'
    txt_f=chr+'_bamview.txt'
    bamview(csv_f,txt_f,bamview_txt_dir)
    # add_pos(chr, txt_dir, txt_add_dir)
    adapter_pic(bam_dir, adapter_dir, chr)

def file_filter(f):
    if f[-4:] in ['.bam']:
        return True
    else:
        return False

if __name__ == '__main__':
    local_time = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print(local_time)
    split_cmd = f"samtools view -H {bam} | cut -f 2 | grep SN | cut -f 2 -d ':' > {txt_dir}/chr.txt"
    with open(log_file, "a") as log:
        os.system(split_cmd)

    chr_list = os.listdir(bam_dir)
    chr_l = list(filter(file_filter, chr_list))
    
    # with tqdm(total=len(chr_l), desc="STAGE 5: SNP-aware adapter removal", bar_format="{l_bar}{bar} |") as pbar:
    #     with ProcessPoolExecutor(max_workers=thread) as executor:  # Changed to ProcessPoolExecutor
    #         futures = [executor.submit(run_adapter, chr) for chr in chr_l]
    #         for future in futures:
    #             future.result()
    #             pbar.update(1)
    if errc==1:
        stage_st=5
    else :
        stage_st=3
    with tqdm(total=len(chr_l), desc=f'STAGE {stage_st}: Adapter detection', bar_format="{l_bar}{bar} |") as pbar:
        with ProcessPoolExecutor(max_workers=thread) as executor:
            futures = [executor.submit(run_adapter, chr) for chr in chr_l]
            # 按任务完成顺序处理
            for future in as_completed(futures):  # 动态监听完成事件
                result = future.result()         # 获取已完成任务的结果
                pbar.update(1)                  # 立即更新进度条              # 立即更新进度条

