import os
import csv
import pysam
import argparse

parser=argparse.ArgumentParser(description='adapter remove.')
parser.add_argument('-outfile',type=str,help='fq with adapter',required=True)
parser.add_argument('-infile',type=str,help='fq',required=True)
parser.add_argument('-csv',type=str,help='csv',required=True)
parser.add_argument('-bam',type=str,help='bam',required=True)
args=parser.parse_args()

outfile = args.outfile
infile = args.infile
csv = args.csv
bam = args.bam

in_file = open(infile,'r') #纠错后的fq文件（含义adapter）

out_file = open(outfile,'w') #去除adapter的fq文件

lines = in_file.readlines()

with open(csv,'r') as f: 
    reader = csv.reader(f)
    rows = list(reader)
del rows[0] #消除第一行

#log_file = open(outdir+"/output.log", "a")
#sys.stdout = log_file
tmp = 0
bamfile = pysam.AlignmentFile(bam, "rb") #bam文件
for bam in bamfile:
    if tmp == len(rows):
        #print('finish',tmp)
        break
    while bam.query_name == rows[tmp][18]:
        row = rows[tmp]
        lenth = round((int(row[5]) - int(row[3]))/3) #adapter数量
        r2 = row[2].split('-')
        r22 = r2[1].split('.')
        #print(r22[0])
        num = int(r22[0]) - 1
        adapter_location = num*445 + int(row[10]) + round(int(row[3])/3) - bam.pos - 1
        count = adapter_location
        true_location = 0 #因为标注框有偏差，所以用bam文件中计算位置
        ci_tmp = 0
        for ci in bam.cigar:
            if ci[0] == 0: #count为比对的位置的计数 如果匹配，则-
                count -= ci[1]
                true_location += ci[1]
            elif ci[0] == 3 or ci[0] == 4:#如果有跳过 locat +
                adapter_location += ci[1]
                true_location += ci[1]
            elif ci[0] == 2 : #如果有删除 
                adapter_location -= ci[1]
                count -= ci[1]
                if(ci[1] <= 5)
                    true_location += ci[1]
            elif ci[0] == 1: #如果有插入
                adapter_location += ci[1]
                true_location += ci[1]
                if(ci[1] <= 5)
                    true_location -= ci[1]
            if  count < 0:
                lenth = bam.cigar[ci_tmp + 1][1]
                if lenth < 30:#代表不是adapter
                    lenth = 0
                break
            ci_tmp += 1
            
        l_r = 0
        #print(true_location,adapter_location,count)
        for i in range(len(lines)):
            if lines[l_r][0] == '@':
                lines[l_r] = lines[l_r].replace('@','')
                line_read_arr = lines[l_r].split()
                lines[l_r] = '@' + lines[l_r]
                if line_read_arr[0] == bam.query_name:
                    l_r += 1
                    str_list = list(lines[l_r])
                    qual_list = list(lines[l_r+2])
                    if bam.flag == 16:
                        """ if len(str_list)  - 2 - read_location < 0 or len(str_list)  - 2 - read_location > len(str_list) - 1:
                            print('error')  """
                        len_str = len(str_list)
                        len_qual = len(qual_list)
                        for j in range(lenth):
                            del str_list[len_str  - 2 - true_location - lenth + 1]  #需要-2 因为最后是/n
                            del qual_list[len_qual  - 2 - true_location - lenth + 1]
                    else:
                        """ if read_location >= len(str_list) or read_location < 0 :
                            print('error',read_location,check_location,len(str_list),read_name,line_ec_arr[0],num)   """
                        for j in range(lenth):
                            del str_list[true_location] 
                            del qual_list[true_location]
                    lines[l_r] = ''.join(str_list)
                    lines[l_r+2] = ''.join(qual_list)
                    break
                else :
                    l_r += 4 #跳到下一条read
            if l_r >= len(lines):
                break

        tmp += 1
        if tmp == len(rows):
            #print('finish',tmp)
            break
        if tmp % 100 == 0:
            #print(tmp)
            now = 1
    
#print(tmp)
out_file.writelines(lines)
#log_file.close()
#sys.stdout = sys.__stdout__ 