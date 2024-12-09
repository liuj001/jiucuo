# -*- coding: utf-8 -*-
import gzip
import os

def get_reads(c_reads_dir, outdir):
    #directory = "/root/autodl-tmp/HG002"
    file_names = list()
    for file in os.listdir(c_reads_dir):
        if os.path.isfile(os.path.join(c_reads_dir, file)):
            file_names.append(os.path.splitext(file))
    print(file_names[0][0]+file_names[0][1])
    print(file_names[0][0])
    print(len(file_names))

    out_file = open(outdir + '/' + 'reads_c.fastq','wt')
    for i in range(len(file_names)):
        in_file = c_reads_dir + '/' + file_names[i][0]+file_names[i][1]
        with gzip.open(in_file, 'rt') as in_f:
            for line in in_f:
                out_file.writelines(line)