import os
import csv
import pysam
import argparse
import gzip
from tqdm import tqdm

parser = argparse.ArgumentParser(description='adapter remove.')
parser.add_argument('-outfile', type=str, help='fq with adapter', required=True)
parser.add_argument('-infile', type=str, help='fq', required=True)
parser.add_argument('-csv', type=str, help='csv', required=True)
parser.add_argument('-bam', type=str, help='bam', required=True)
parser.add_argument('-error_correction', type=int, help='had corrected error', required=True)
args = parser.parse_args()

outfile = args.outfile
infile = args.infile
csv_f = args.csv
bam = args.bam
isonly = 0 # 0: remove adapter after error correction 1: remove adapter only
if args.error_correction == 0:
    isonly = 1

in_file = open(infile, 'r')  # 纠错后的fq文件（含义adapter）
out_file = open(outfile, 'w')  # 去除adapter的fq文件

lines = in_file.readlines()

with open(csv_f, 'r') as f: 
    reader = csv.reader(f)
    rows = list(reader)


header = rows[0]
data = rows[1:]

# 按第一列升序排序（不区分大小写）
sorted_data = sorted(data, key=lambda x: x[0].lower())

# 合并标题与排序后的数据
rows = [header] + sorted_data
del rows[0]  # 消除第一行

tmp = 0
bamfile = pysam.AlignmentFile(bam, "rb")  # bam文件

# 初始化进度条
if isonly == 1:
    stage_st=4
else :
    stage_st=6

with tqdm(desc="STAGE {stage_st}: Adapter removal", bar_format="{l_bar}{bar} |") as pbar:
    pbar.total = 0  # 初始化总量为 0
    for bam in bamfile:
        pbar.total += 1  # 动态增加总量
        pbar.update(1)  # 更新进度条

        if tmp == len(rows):
            break
        while bam.query_name == rows[tmp][16]:
            row = rows[tmp]
            lenth = round((int(row[3]) - int(row[1])) / 3)  # adapter数量
            r2 = row[0].split('-')
            r22 = r2[1].split('.')
            num = int(r22[0]) - 1
            adapter_location = num * 445 + int(row[8]) + round(int(row[1]) / 3) - bam.pos - 1
            count = adapter_location
            true_location = 0  # 因为标注框有偏差，所以用bam文件中计算位置
            ci_tmp = 0
            for ci in bam.cigar:
                if ci[0] == 0:  # count为比对的位置的计数 如果匹配，则-
                    count -= ci[1]
                    true_location += ci[1]
                elif ci[0] == 3 or ci[0] == 4:  # 如果有跳过 locat +
                    adapter_location += ci[1]
                    true_location += ci[1]
                elif ci[0] == 2:  # 如果有删除 
                    adapter_location -= ci[1]
                    count -= ci[1]
                    if ci[1] <= 5 and isonly == 0:
                        true_location += ci[1]
                elif ci[0] == 1:  # 如果有插入
                    adapter_location += ci[1]
                    true_location += ci[1]
                    if ci[1] <= 5 and isonly == 0:
                        true_location -= ci[1]
                if count < 0:
                    lenth = bam.cigar[ci_tmp + 1][1]
                    if lenth < 30:  # 代表不是adapter
                        lenth = 0
                    break
                ci_tmp += 1

            l_r = 0
            if lenth != 0:
                for i in range(len(lines)):
                    if lines[l_r][0] == '@':
                        lines[l_r] = lines[l_r].replace('@', '')
                        line_read_arr = lines[l_r].split()
                        lines[l_r] = '@' + lines[l_r]
                        if line_read_arr[0] == bam.query_name:
                            l_r += 1
                            str_list = list(lines[l_r])
                            qual_list = list(lines[l_r + 2])
                            if bam.flag == 16:
                                len_str = len(str_list)
                                len_qual = len(qual_list)
                                for j in range(lenth):
                                    if str_list[len_str - 2 - true_location - lenth + 1] != '\n':
                                        del str_list[len_str - 2 - true_location - lenth + 1]  # 需要-2 因为最后是/n
                                        del qual_list[len_qual - 2 - true_location - lenth + 1]
                            else:
                                for j in range(lenth):
                                    if str_list[true_location] != '\n':
                                        del str_list[true_location]
                                        del qual_list[true_location]
                            del str_list[-1]
                            del qual_list[-1]
                            seq_l = ''.join(str_list[:true_location]) + '\n'
                            qual_l = ''.join(qual_list[:true_location]) + '\n'
                            seq_r = ''.join(str_list[true_location:]) + '\n'
                            qual_r = ''.join(qual_list[true_location:]) + '\n'
                            if true_location < len(str_list)/10:
                                lines[l_r] = seq_r
                                lines[l_r + 2] = qual_r
                            elif true_location > len(str_list)*9/10:
                                lines[l_r] = seq_l
                                lines[l_r + 2] = qual_l
                            else:
                                lines[l_r-1] = bam.query_name + '_1\n'
                                lines[l_r] = seq_l
                                lines[l_r + 2] = qual_l
                                lines[l_r+3] = bam.query_name + '_2\n'
                                lines[l_r+4] = seq_r
                                lines[l_r+5] = '+\n'
                                lines[l_r+6] = qual_r

                            # lines[l_r] = ''.join(str_list)
                            # lines[l_r + 2] = ''.join(qual_list)
                            break
                        else:
                            l_r += 4  # 跳到下一条read
                    if l_r >= len(lines):
                        break

            tmp += 1
            if tmp == len(rows):
                break

out_file.writelines(lines)
