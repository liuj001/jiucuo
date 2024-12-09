import pandas as pd
from sklearn.cluster import DBSCAN

import numpy as np
import regex as re
from datetime import datetime
import os
import csv
import pysam
import argparse
from collections import Counter
# need detection_result.csv


def find_substring_indices(long_string, string_array,window_size):
    substring_indices = {}
    substring_length = window_size
    for i in range(len(long_string) - substring_length + 1):
        substring = long_string[i:i + substring_length]
        if substring in string_array:
            if substring not in substring_indices:
                substring_indices[substring] = []
            substring_indices[substring].append(i)
    return substring_indices


def process_dict(input_dict):
    # 将字典的值转换为数组列表
    array_list = list(input_dict.values())

    # 初始化结果数组
    result = []

    for i in range(len(array_list)):
        current_array = array_list[i]
        if i == 0:
            # 第一个数组直接选择第一个数字
            result.append(current_array[0])
        else:
            previous = result[-1]
            found = False
            for num in current_array:
                if num == previous + 1:
                    result.append(num)
                    found = True
                    break
            if not found:
                for num in current_array:
                    if num > previous:
                        result.append(num)
                        break

    return result


def calculate_result(window_size, arr):
    result = 0
    n = len(arr)
    for i in range(n):
        if i==0:
            result+=window_size
        elif arr[i]-arr[i-1] >= window_size:
            result+=window_size
        elif arr[i]-arr[i-1]<window_size:
            x=arr[i]-arr[i-1]
            result+=x
    return result

def union_func(long_string, string_array,window_size,len):
    indices=find_substring_indices(long_string, string_array,window_size)
    arr=process_dict(indices)
    calculated_result=calculate_result(window_size, arr)
    return round(calculated_result/len,2)



def contains_adapter_check(x, window_size, thresholds,new_array1,new_array2,new_array3,new_array4):

    for new_array, threshold in thresholds.items():
        count=0
        if len(x)==35:
            
            if(new_array=='new_array2'):
                count=union_func(x,new_array2,window_size,len=35)
            if(new_array=='new_array4'):
                count=union_func(x,new_array4,window_size,len=35)
        if len(x)==45:
            if(new_array=='new_array1'):
                count=union_func(x,new_array1,window_size,len=45)
            if(new_array=='new_array3'):
                count=union_func(x,new_array3,window_size,len=45)
      
        if count >= threshold:
            if count>=1:
                count=1
            print("count:"+str(count)+" x:"+str(x)+" threshold"+str(threshold))
            
       
            return True

    return False

def base_all_detections1(openfile,bam_filename):
    out_file = open(f'{openfile}_pre_adapter_out_all_detections.csv','w')
    # 使用 pandas 读取 CSV 文件为 DataFrame  
    # df1 = pd.read_csv(f'{openfile}_p.csv')  
    df1 = pd.read_csv(f'{openfile}_p.csv')  

    df1.sort_values(by='trueY1', ascending=True,inplace=True)  
 
  
    df1.to_csv(f'{openfile}_p_sorted_file.csv', index=True)
    with open(f'{openfile}_p_sorted_file.csv','r') as in_file:
        reader = csv.reader(in_file)
        rows = list(reader)
      
    writer = csv.writer(out_file)
    rows[0].append('read_seq')
    writer.writerow(rows[0])
    tmp = 1
    bamfile = pysam.AlignmentFile(bam_filename, "rb")

    for bam in bamfile:
        """ if bam.pos == 60586:
            print(bam.cigar)
            break """
        while bam.query_name == rows[tmp][18]:
            
            row = rows[tmp]
            start = int(row[10])
            r2 = row[2].split('-')
            r22 = r2[2].split('.')
           
            num = int(r22[0]) - 1
            truex1 = start + num*445 + int( (int(row[3]) -1) / 3 )
            truex2 = start + num*445 + int( (int(row[5]) -1) / 3 )
            count = 0
            locat = 0
            for ci in bam.cigar:
                if ci[0] == 4 or ci[0] == 1:
                    count += ci[1]
                elif ci[0] == 2:
                    count -= ci[1]
                    locat += ci[1]
                elif ci[0] == 0:
                    locat += ci[1]
                if locat >= truex1 - bam.pos:
                    break
            if truex1 - bam.pos <= 1:#头
                count -= 45
            elif truex2-bam.pos >= locat + 30:#尾
                count += 1
            truex1 += count
            truex2 += count
            if truex1 - bam.pos < 0:
                add = bam.pos-truex1
                truex1 += add
                truex2 += add
            if truex1-bam.pos-3>=0:
                left=truex1-bam.pos
            else:
                left=truex1-bam.pos
            sequence = bam.query_sequence[left:truex2-bam.pos+3]
            
            row.append(sequence)
            writer = csv.writer(out_file)
            writer.writerow(row)
            tmp += 1
            if tmp == len(rows):
            
                break
        if tmp == len(rows):
            break

def ATGCFiltDetections1(openfile,similarity,k_size):
    # 读取CSV文件
    df = pd.read_csv(f'{openfile}_pre_adapter_out_all_detections.csv')

    # 定义四种adapter序列
    adapter_sequences = ['ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT', 'AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA',
                         'ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT', 'TCCTCCTCCTCCGTTAATTTTTTTTTTTTTTTTTT']
    # adapter_sequences_houxuan= [, , , ]
    new_array1 = []
    new_array2 = []
    new_array3 = []
    new_array4 = []
    window=k_size
    sequence0 = adapter_sequences[0]
    for i in range(len(sequence0) - window+1):
        new_array1.append(sequence0[i:i + window])

    sequence1 = adapter_sequences[1]
    for i in range(len(sequence1) - window+1):
        new_array2.append(sequence1[i:i + window])

    sequence2 = adapter_sequences[2]
    for i in range(len(sequence2) - window+1):
        new_array3.append(sequence2[i:i + window])

    sequence3 = adapter_sequences[3]
    for i in range(len(sequence3) - window +1):
        new_array4.append(sequence3[i:i + window])

    df['read_seq'] = df['read_seq'].astype(str)
    



    thresholds = {'new_array1':similarity, 'new_array2':similarity, 'new_array3':similarity, 'new_array4':similarity}
    window_size = 8


    size=45
    df['contains_adapter'] = df['read_seq'].apply(
    lambda seq: any(
        contains_adapter_check(sub_seq, window_size, thresholds, new_array1, new_array2, new_array3, new_array4)
        for sub_seq in [seq[i:i+size] for i in range(len(seq) - size + 1)]
    )
    )
   
    size1=35
    df['contains_adapter1'] = df['read_seq'].apply(
    lambda seq: any(
        contains_adapter_check(sub_seq, window_size, thresholds, new_array1, new_array2, new_array3, new_array4)
        for sub_seq in [seq[i:i+size1] for i in range(len(seq) - size1 + 1)]
    )
    )

    adapter_rows = df[df['contains_adapter'] == True]
    adapter_rows1 = df[df['contains_adapter'] == True]

    file_name = f'all_detections_to_pre_adapter.csv'
    adapter_rows_columns_to_save=adapter_rows.columns[2:]
    data_to_save=adapter_rows[adapter_rows_columns_to_save]


    adapter_rows1_columns_to_save=adapter_rows1.columns[2:]
    data_to_save1=adapter_rows1[adapter_rows1_columns_to_save]
   

    combined_data = pd.concat([data_to_save1, data_to_save])
    combined_data = combined_data.drop_duplicates(subset=['image_id', 'x1'])
    
    combined_data.to_csv(f'{openfile}_{file_name}',index=False)
    openfile=openfile[8:-4]
    if not combined_data.empty:
        combined_data.to_csv(f'../result/{openfile}_{file_name}',index=False)
    print(f'../result/{openfile}_{file_name}')
    
    
    
 

    return adapter_rows


