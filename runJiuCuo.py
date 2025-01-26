# -*- coding: utf-8 -*-
import os
import argparse
import pysam
from time import *
from concurrent.futures import ThreadPoolExecutor
from threading import Thread
import subprocess
import sys
import logging

from tqdm import tqdm
from picture.candidate_make import candidate_make
from picture.deepvariant_speed_final_ex1 import snp_pic
from picture.bamview_new import bamview
from picture.add_pos import add_pos
from picture.Ex_k4_adapter_ver2 import adapter_pic
from inception.predict_a import find_snp
from correction.find_error import find_error
from correction.read_locate_indel_v5 import correct

import glob

parser=argparse.ArgumentParser(description='runJiuCuo.py integrates the processes of error correction.')
parser.add_argument('-contigs',type=str,help='Preassembled primary contigs from the reads in FASTA format',required=True)
parser.add_argument('-reads',type=str,help='Raw HiFi reads in FASTQ format',required=True)
parser.add_argument('-output',type=str,help='Output directory',required=True)
#parser.add_argument('-min_bases',type=int,help='Minumum number of mismatched bases required to generate an error candidate image',default=1)
#parser.add_argument('-min_reads',type=int,help='Minimum number of reads required to generate the error candidate image',default=3)
parser.add_argument('-threads',type=int,help='Number of threads during correction',default=8)
parser.add_argument('-allocated_reads',type=int,help='Maximum number of reads allocated to each thread',default=10000)
parser.add_argument('-adapter_removal',type=int,help='Adapter removal from the reads.',default=0)
args=parser.parse_args()

ref = args.contigs
reads_path = args.reads
outdir = args.output
#min_bases = args.min_bases
#min_reads = args.min_reads
thread = args.threads
allocated_reads = args.allocated_reads
adapter_removal = args.adapter_removal
bam = outdir+'/raw.bam'
bam_dir = outdir+'/bam'
txt_dir = outdir+'/txt'
txt_add_dir = outdir+'/txt_add'
bcf_txt_dir = outdir+'/bcf_txt'
snp_c_dir = outdir+'/snp_candidate'
adapter_dir = outdir+'/adapter'
adapter_out_dir = outdir+'/adapter_out'
bamview_txt_dir = outdir+'/bamview_txt'
bamview_csv_dir = outdir+'/bamview_csv'
snp_dir = outdir+'/snp'
ec_dir = outdir+'/ec_txt'
c_reads_dir = outdir+'/c_reads'
# 日志文件
log_file = outdir+'/runJiuCuo.log'

'''make outdir'''
if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(bam_dir):
    os.makedirs(bam_dir)
if not os.path.exists(txt_dir):
    os.makedirs(txt_dir)
if not os.path.exists(txt_add_dir):
    os.makedirs(txt_add_dir)
if not os.path.exists(bcf_txt_dir):
    os.makedirs(bcf_txt_dir)
if not os.path.exists(snp_c_dir):
    os.makedirs(snp_c_dir)
if not os.path.exists(adapter_dir):
    os.makedirs(adapter_dir)
if not os.path.exists(adapter_out_dir):
    os.makedirs(adapter_out_dir)
if not os.path.exists(bamview_txt_dir):
    os.makedirs(bamview_txt_dir)
if not os.path.exists(bamview_csv_dir):
    os.makedirs(bamview_csv_dir)
if not os.path.exists(snp_dir):
    os.makedirs(snp_dir)
if not os.path.exists(ec_dir):
    os.makedirs(ec_dir)
if not os.path.exists(c_reads_dir):
    os.makedirs(c_reads_dir)

def re_refname(chr):
    inbam = os.path.join(bam_dir, f'{chr}.bam')
    in_bam = pysam.AlignmentFile(inbam,'rb')
    header = {'SQ': [{'LN': 0, 'SN': chr}]}
    outbam = os.path.join(bam_dir, f'{chr}_reheader.bam')
    out_bam = pysam.AlignmentFile(outbam,'wb',header = header)
    for bam in in_bam:
        a = pysam.AlignedSegment()
        a = bam
        a.reference_id = 0
        out_bam.write(a)
    in_bam.close()
    out_bam.close()
    os.remove(inbam)
    os.rename(outbam, inbam)

def error_correct(chr):
    start = time()
    chr = chr.rstrip('.bam')
    reads_dir_path = c_reads_dir+'/'+chr+'.ec.fastq.gz'
    if not os.path.exists(reads_dir_path):
        if not os.path.exists(bam_dir+'/'+chr+'.bam.bai'):
            bai_cmd = "samtools index %s/%s.bam %s/%s.bam.bai"%(bam_dir, chr, bam_dir, chr)
            with open(log_file, "a") as log:
                bai_p=subprocess.Popen(bai_cmd,shell=True, stdout=log, stderr=log)
                bai_code=bai_p.wait()  #等待子进程结束，并返回状态码;
            #samtools index -@ 48 -b /root/autodl-tmp/HG002/old/HG002.raw.sort.filt.bam
            #os.system("samtools index {out}/{chr}.bam {out}/{chr}.bam.bai".format(out=bam_dir,chr=chr))
        if not os.path.exists(txt_dir+'/'+chr+'.txt'):
            txt_cmd = "samtools mpileup -d 100000 -Q 0 -f %s %s/%s.bam > %s/%s.txt"%(ref, bam_dir, chr, txt_dir, chr)
            with open(log_file, "a") as log:
                txt_p=subprocess.Popen(txt_cmd,shell=True, stdout=log, stderr=log)
                txt_code=txt_p.wait()  #等待子进程结束，并返回状态码;
        
        snp_pic_dir = snp_c_dir+'/'+chr
        txt_file = txt_dir + '/' + chr + '.txt'
        sz = os.path.getsize(txt_file)
        if not sz:
            end = time()
            #print(chr,"done!", f'time: {end - start:.3f}s.=========================')
        else:  
            if adapter_removal==1:
                bamview_txt_cmd = "samtools view %s/%s.bam > %s/%s_bamview.txt"%(bam_dir, chr, bamview_txt_dir, chr)
                with open(log_file, "a") as log:
                    bamview_txt_p=subprocess.Popen(bamview_txt_cmd,shell=True, stdout=log, stderr=log)
                    bamview_txt_code=bamview_txt_p.wait()  #等待子进程结束，并返回状态码;
                csv_f=bamview_csv_dir+'/'+chr+'.csv'
                txt_f=chr+'_bamview.txt'
                bamview(csv_f,txt_f,bamview_txt_dir)
                add_pos(chr, txt_dir, txt_add_dir)
                adapter_pic(bam_dir, txt_add_dir, adapter_dir, chr)

            #print(chr,"make candidate =========================")  
            candidate_make(chr, txt_dir, bcf_txt_dir)
            bcf_file = bcf_txt_dir + '/' + chr + '.bcf.txt'
            szbcf = os.path.getsize(bcf_file)
            if not szbcf:
                #print(chr,"has no snp candidate")
                nosnp=1
            else:
                #print(chr,"pic =========================")
                if not os.path.exists(snp_pic_dir):
                    os.makedirs(snp_pic_dir)
                snp_pic(bam_dir, txt_dir, bcf_txt_dir, snp_pic_dir, chr)
                #print(chr,"call snp =========================")
                find_snp(snp_pic_dir, snp_dir, chr)
            
            #print(chr,"find error =========================")
            find_error(ec_dir, snp_dir, txt_dir, chr)
            #print(chr,"correct error =========================")
            correct(bam_dir, txt_dir, ec_dir, c_reads_dir, chr)
            end = time()
            re_refname(chr)
            #print(chr,"done!", f'time: {end - start:.3f}s.=========================')
    
    return chr

def spilt(chr):
    chr = chr.rstrip('\n')
    start = time()
    if not os.path.exists(bam_dir+'/'+chr+'_1.bam'):
        if not os.path.exists(bam_dir+'/'+chr+'.bam'):
            bam_cmd = "samtools view -@ 12 -b %s/filt.bam %s > %s/%s.bam"%(outdir, chr, bam_dir, chr)
            with open(log_file, "a") as log:
                bam_p=subprocess.Popen(bam_cmd,shell=True, stdout=log, stderr=log)
                bam_code=bam_p.wait()  #等待子进程结束，并返回状态码;
            #samtools view -@ 12 -b test.bam chr1 > chr1.bam

        split_size = allocated_reads
        # 打开输入 BAM 文件
        input_bam = bam_dir+'/'+chr+'.bam'
        bamfile = pysam.AlignmentFile(input_bam, "rb")
        
        # 获取 BAM 文件的头信息，用于创建新的 BAM 文件
        header = bamfile.header
        
        # 初始化变量
        counter = 0
        split_counter = 0
        current_bam_writer = None

        for read in bamfile:
            # 每 split_size 个位置创建一个新的 BAM 文件
            if counter % split_size == 0:
                if current_bam_writer:
                    current_bam_writer.close()  # 关闭之前的 BAM 文件
                split_counter += 1
                out_bam = f"{chr}_{split_counter}.bam"
                output_bam = bam_dir+'/'+ out_bam
                current_bam_writer = pysam.AlignmentFile(output_bam, "wb", header=header)

            # 将当前读取的 read 写入当前的 BAM 文件
            current_bam_writer.write(read)
            
            # 增加计数器
            counter += 1

        # 关闭最后一个 BAM 文件
        if current_bam_writer:
            current_bam_writer.close()

        bamfile.close()
        #print(f"Total {split_counter} BAM files created.")
    
        os.remove(input_bam)
    end = time()
    #print(chr,"spilt done!", f'time: {end - start:.3f}s.=========================') 
    
    return chr

# def main():
#     if not os.path.exists(ref+'.fai'):
#         #print("faidx")
#         os.system("samtools faidx {ref}".format(ref=ref))
#         #print("Done!")
#     #samtools faidx /root/autodl-tmp/HG002/hifiasm/m64012_190920_173625.raw.bp.p_ctg.fa
#     if not os.path.exists(outdir+'/sort.bam'):
#         #print("Sort bam")
#         sort_cmd = "samtools sort -@ 48 %s > %s/sort.bam"%(bam,outdir)
#         sort_p=subprocess.Popen(sort_cmd,shell=True)
#         sort_code=sort_p.wait()  #等待子进程结束，并返回状态码；
#         #print("Done!")
#     #samtools sort -@ 48 -o /root/autodl-tmp/test/A.tha_hifiasm/CRR302668.asm.sort.bam /root/autodl-tmp/test/A.tha_hifiasm/CRR302668.asm.bam
#     if not os.path.exists(outdir+'/filt.bam'):
#         #print("Filt bam")
#         filt_cmd = "samtools view -@ 48 -b -F 4 -F 256 -F 2048 %s/sort.bam > %s/filt.bam"%(outdir,outdir)
#         filt_p=subprocess.Popen(filt_cmd,shell=True)
#         filt_code=filt_p.wait()  #等待子进程结束，并返回状态码；
#         #print("Done!")
#     #samtools view -@ 48 -b -F 4 -F 256 -F 2048 /root/autodl-tmp/HG002/m64012_190920_173625.raw.sort.bam > /root/autodl-tmp/HG002/m64012_190920_173625.raw.sort.filt.bam
#     if not os.path.exists(outdir+'/filt.bam.bai'):
#         #print("Index bam")
#         bai_cmd = "samtools index %s/filt.bam %s/filt.bam.bai"%(outdir,outdir)
#         bai_p=subprocess.Popen(bai_cmd,shell=True)
#         bai_code=bai_p.wait()  #等待子进程结束，并返回状态码;
#         #print("Done!")
#     #samtools index -@ 48 -b /root/autodl-tmp/HG002/old/HG002.raw.sort.filt.bam
#     split_cmd = "samtools view -H %s | cut -f 2 | grep SN | cut -f 2 -d \":\" > %s/chr.txt"%(bam,txt_dir)
#     os.system(split_cmd) #得到该bam文件的所有染色体号

def main():
    if not os.path.exists(ref + '.fai'):
        faidx_cmd = f"samtools faidx {ref}"
        with open(log_file, "a") as log:
            faidx_p = subprocess.Popen(faidx_cmd, shell=True, stdout=log, stderr=log)
            faidx_p.wait()

    if not os.path.exists(outdir + '/sort.bam'):
        sort_cmd = f"samtools sort -@ 48 {bam} > {outdir}/sort.bam"
        with open(log_file, "a") as log:
            sort_p = subprocess.Popen(sort_cmd, shell=True, stdout=log, stderr=log)
            sort_p.wait()

    if not os.path.exists(outdir + '/filt.bam'):
        filt_cmd = f"samtools view -@ 48 -b -F 4 -F 256 -F 2048 {outdir}/sort.bam > {outdir}/filt.bam"
        with open(log_file, "a") as log:
            filt_p = subprocess.Popen(filt_cmd, shell=True, stdout=log, stderr=log)
            filt_p.wait()

    if not os.path.exists(outdir + '/filt.bam.bai'):
        index_cmd = f"samtools index {outdir}/filt.bam {outdir}/filt.bam.bai"
        with open(log_file, "a") as log:
            index_p = subprocess.Popen(index_cmd, shell=True, stdout=log, stderr=log)
            index_p.wait()

    split_cmd = f"samtools view -H {bam} | cut -f 2 | grep SN | cut -f 2 -d ':' > {txt_dir}/chr.txt"
    with open(log_file, "a") as log:
        os.system(split_cmd)

def file_filter(f):
    if f[-4:] in ['.bam']:
        return True
    else:
        return False

if __name__ == '__main__':
    start = time()
    main()
    end_bam = time()
    #print(f'==============================Filt bam: {end_bam - start:.3f}s.==============================')
    
    #start_spilt = time()
    #with open(txt_dir+'/chr.txt','r') as chr_file:
    #    chr_n = chr_file.readlines()
    #with ThreadPoolExecutor(max_workers=thread) as executor:
    #    res = executor.map(spilt, chr_n)

    #end_spilt = time()
    #print(f'==============================Split bam: {end_spilt - start_spilt:.3f}s.==============================')
    start_spilt = time()
    with open(txt_dir + '/chr.txt', 'r') as chr_file:
        chr_n = chr_file.readlines()

    # 创建进度条
    with tqdm(total=len(chr_n), desc="Splitting BAM files") as pbar:
        with ThreadPoolExecutor(max_workers=thread) as executor:
            # 使用回调函数更新进度条
            futures = [executor.submit(spilt, chr) for chr in chr_n]
            for future in futures:
                future.result()
                pbar.update(1)

    end_spilt = time()

    # start_error = time()
    # chr_list= os.listdir(bam_dir)
    # chr_l = list(filter(file_filter, chr_list))
    # with ThreadPoolExecutor(max_workers=thread) as executor:
    #     res = executor.map(error_correct, chr_l)
    # end_error = time()
    #print(f'==============================Correct error: {end_error - start_error:.3f}s.==============================')    
    start_error = time()
    chr_list = os.listdir(bam_dir)
    chr_l = list(filter(file_filter, chr_list))

    # 创建进度条
    with tqdm(total=len(chr_l), desc="Correcting Errors") as pbar:
        with ThreadPoolExecutor(max_workers=thread) as executor:
            futures = [executor.submit(error_correct, chr) for chr in chr_l]
            for future in futures:
                future.result()
                pbar.update(1)

    end_error = time()
     
    start_read = time()
    # 显示拼接进度
    print("Processing reads...")
    with tqdm(total=3, desc="Merging and Processing Reads") as pbar:
        read_cmd = "cat %s/*.gz > %s/reads_c.fastq.gz" % (c_reads_dir, outdir)
        with open(log_file, "a") as log:
            read_p = subprocess.Popen(read_cmd, shell=True, stdout=log, stderr=log)
            read_code = read_p.wait()
        pbar.update(1)

        com_cmd = "seqkit common %s %s/reads_c.fastq.gz | seqkit seq -n -i -o %s/common.txt" % (reads_path, outdir, outdir)
        with open(log_file, "a") as log:
            com_p = subprocess.Popen(com_cmd, shell=True, stdout=log, stderr=log)
            com_code = com_p.wait()
        pbar.update(1)

        rem_cmd = "seqkit grep -v -f %s/common.txt %s -o %s/s.fastq.gz" % (outdir, reads_path, outdir)
        with open(log_file, "a") as log:
            rem_p = subprocess.Popen(rem_cmd, shell=True, stdout=log, stderr=log)
            rem_code = rem_p.wait()
        pbar.update(1)

    t = outdir + "/common.txt"
    os.remove(t)

    end_cmd = "cat %s/s.fastq.gz %s/reads_c.fastq.gz > %s/base_correction.fastq.gz" % (outdir, outdir, outdir)
    with open(log_file, "a") as log:
        end_p = subprocess.Popen(end_cmd, shell=True, stdout=log, stderr=log)
        end_code = end_p.wait()

    c = outdir + "/reads_c.fastq.gz"
    os.remove(c)
    s = outdir + "/s.fastq.gz"
    os.remove(s)

    end = time()


    # start_read = time()
    # #print("Get reads ====================")
    # read_cmd = "cat %s/*.gz > %s/reads_c.fastq.gz"%(c_reads_dir, outdir) 
    # read_p=subprocess.Popen(read_cmd,shell=True)
    # read_code=read_p.wait()

    # #print("Find common====================")
    # com_cmd = "seqkit common %s %s/reads_c.fastq.gz | seqkit seq -n -i -o %s/common.txt"%(reads_path, outdir, outdir)
    # #seqkit common file*.fa | seqkit seq -n -i -o common.txt
    # com_p=subprocess.Popen(com_cmd,shell=True)
    # com_code=com_p.wait()

    # #print("Remove====================")
    # rem_cmd = "seqkit grep -v -f %s/common.txt %s -o %s/s.fastq.gz"%(outdir, reads_path, outdir)
    # #seqkit grep -v -f common.txt a.fa -o common.txt
    # rem_p=subprocess.Popen(rem_cmd,shell=True)
    # rem_code=rem_p.wait()

    # t = outdir+"/common.txt"
    # os.remove(t)

    # #print("Cat====================")
    # end_cmd = "cat %s/s.fastq.gz %s/reads_c.fastq.gz > %s/correction.fastq.gz"%(outdir, outdir, outdir)
    # #cat Sample_test_1.R1.fastq.gz Sample_test_2.R2.fastq.gz > test2.fastq.gz
    # end_p=subprocess.Popen(end_cmd,shell=True)
    # end_code=end_p.wait()

    # c = outdir+"/reads_c.fastq.gz"
    # os.remove(c)
    # s = outdir+"/s.fastq.gz"
    # os.remove(s)

    # end = time()
    # #print(f'==============================Total time: {end - start_read:.3f}s.==============================')
    # #print(f'==============================Total time: {end - start:.3f}s.==============================')