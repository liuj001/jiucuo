#snp生成图片加快版本
import numpy as np
import pysam
from numpy import * 
from PIL import Image
from time import *
import json
import torch
from itertools import islice

import os

#path_bam = 'chr20_d70.asm.bp.p_ctg.bam' #bam
#path_txt = 'chr20_d70.asm.bp.p_ctg.txt' #txt
#path_pic = 'pic/input/snp'  #图片路径
#path_bcf = 'chr20_d70_bcf.txt'  #候选文件路径
#path_gt = '../SNPTools/benchmark/20.vcf'

def extract_before_underscore(strings):
    return ['_'.join(s.split('_')[:-1]) for s in strings]

def snp_pic(bam_dir, txt_dir, bcf_txt_dir, snp_pic, chr):
    #print("snp_pic")
    path_bam = bam_dir + '/' + chr + '.bam'
    path_txt = txt_dir + '/' + chr + '.txt'
    path_bcf = bcf_txt_dir + '/' + chr + '.bcf.txt'
    path_pic = snp_pic
    begin_time = time()
    ## 获取reads的名称
    #bamfile = pysam.AlignmentFile(path_bam, "rb")
    refs_ = []
    refs_.append(chr)

    refs = extract_before_underscore(refs_)
    #bamfile.close()
    # print('Refrence Done!')
    ## 读取比对文件
    #begin_txt = time()
    #txtfile = open(path_txt,"r")

    #lines = txtfile.readlines() 

    #txtfile.close()
    #print('1')
    bcffile = open(path_bcf,"r") #打开bcf文件
    bcf_lines = bcffile.readlines()
    bcffile.close()
    #print('2')
    #gtfile = open(path_gt,'r')
    #gt_lines = gtfile.readlines()
    """ gt_lines = list()
    with open(path_gt, 'r') as gtfile:
        for gt_line in gtfile:
            gt_line_arr = gt_line.split()
            if gt_line_arr[0] == 'chr20' :
                gt_lines.append(gt_line)
            elif gt_line_arr[0] == 'chr21' :
                break 
    #gtfile.close()
    print('txt Done!')
    end_txt = time()
    run_txt = end_txt-begin_txt
    print ('txt读取时间：',run_txt) """

    ref_n = 0
    bcf_n = 1
    ## 以每一条参考read为一组
    run_time = time()
    for ref in refs:
        lines = list()
        lines_len = 0
        begin_txt = time()
        with open(path_txt, 'r') as f:
            for line in islice(f,0,None):
                line_arr = line.split()
                if line_arr[0] == ref:
                    lines.append(line)
                    lines_len += 1
                elif lines_len != 0:
                    break
        #print('txt Done!')
        end_txt = time()
        run_txt = end_txt-begin_txt
        #print ('txt读取时间：',run_txt)
        f.close()
        if ref == '':
            break
        else :
            #print('1111')
            #print(lines_len)
            num_snp = 0
            num_nosnp = 0
            num_skip = 0
            num_skipp = 0
            begin_r = time()
            #print(ref)
            c = 0
            #c1 = 0
            img_n = 1
            ref_n += 1
            count = 0
            #count_bcf = 0
            #count_bcf1 = 0
            #head = 0
            jump = 0
            #jump_deep = 0 #仿照deepvariant每张图第一行都有信息
            #flag_prior = 0
            #flag_break = 0
            #img_a = array([[0 for x in range(2000)] for y in range(1200)],dtype=uint8) #actg数组
            #img_b = array([[0 for x in range(2000)] for y in range(1200)],dtype=uint8) #测序质量
            #img_c = array([[0 for x in range(2000)] for y in range(1200)],dtype=uint8) #对错情况
            #img_a1 = array([[0 for x in range(60)] for y in range(1200)],dtype=uint8)
            #img_b1 = array([[0 for x in range(60)] for y in range(1200)],dtype=uint8)
            #img_c1 = array([[0 for x in range(60)] for y in range(1200)],dtype=uint8)
            #img_mapping_quality = array([[0 for x in range(31)] for y in range(1200)],dtype=uint8)
            #img_strand = array([[0 for x in range(31)] for y in range(1200)],dtype=uint8)
            img_read_base = array([[0 for x in range(31)] for y in range(32)],dtype=uint8)
            img_base_quality = array([[0 for x in range(31)] for y in range(32)],dtype=uint8)
            img_mq = array([[0 for x in range(31)] for y in range(32)],dtype=uint8)
            img_st = array([[0 for x in range(31)] for y in range(32)],dtype=uint8)
            img_conbine = array([[0 for x in range(186)] for y in range(32)],dtype=uint8)
            img_bdfr = array([[0 for x in range(31)] for y in range(32)],dtype=uint8)
            img_rsv = array([[7 for x in range(31)] for y in range(32)],dtype=uint8)
            pos = list()
            tmp_pos = 0
            #记录所有snp候选位置为pos
            for bcf_line in bcf_lines:
                if bcf_line[0] != '#':
                    bcf_arr = bcf_line.split()
                    if bcf_arr[0] == ref:
                        pos.append(bcf_arr[1])
            #print(len(pos))
            if len(pos) == 0:
                pass
            else:
                """ gt = list()
                tmp_gt = 0
                #记录所有gt位置
                for gt_line in gt_lines:
                    if gt_line[0] != '#':
                        gt_arr = gt_line.split()
                        if gt_arr[0] == ref:
                            gt.append(gt_arr[1]) """
                #记录bam文件中的flag,mapping qualty到数组里
                bam_array = list()
                bamfile = pysam.AlignmentFile(path_bam, "rb")
                time1 = time()
                a = 0
                for bam in bamfile:
                    if bam.reference_name != ref and len(bam_array) != 0:
                        break
                    if (bam.flag != 4 and bam.flag != 256 and bam.flag != 2048) and bam.reference_name == ref :
                        bam_array.append([bam.flag,bam.mapping_quality])
                bamfile.close()
                #print('bam done',len(bam_array))
                time2 = time()
                #print(time2 - time1)

                for line in lines:
                    count += 1
                    line_arr = line.split()
                    if line_arr[0] != ref:
                        #print('error')
                        break
                    basess = list(line_arr[4])
                    while int(line_arr[1])+15 > int(pos[tmp_pos]):
                        tmp_pos += 1
                        num_skip += 1
                        if tmp_pos == len(pos):
                            break
                    if int(line_arr[1])+15 == int(pos[tmp_pos])  : #当遍历到snp的区域，进行数组赋值
                        flag_pos = pos[tmp_pos]
                        img_read_base = array([[0 for x in range(31)] for y in range(32)],dtype=uint8)
                        img_base_quality = array([[0 for x in range(31)] for y in range(32)],dtype=uint8)
                        img_mq = array([[0 for x in range(31)] for y in range(32)],dtype=uint8)
                        img_st = array([[0 for x in range(31)] for y in range(32)],dtype=uint8)
                        img_conbine = array([[0 for x in range(186)] for y in range(32)],dtype=uint8)
                        img_bdfr = array([[0 for x in range(31)] for y in range(32)],dtype=uint8)
                        img_rsv = array([[7 for x in range(31)] for y in range(32)],dtype=uint8)
                        b = time()
                        for line_31 in lines[count-1:count+30]:
                            line_31_arr = line_31.split()
                            bases = list(line_31_arr[4])
                            qual = list(line_31_arr[5])
                            b_len = len(bases)
                            q_len = len(qual)
                            # print(b_len)
                            bases_a = array([0 for x in range(b_len)])
                            qual_a = array([0 for x in range(b_len)])
                            bdfr_a = array([0 for x in range(b_len)])
                            l_b = 0 
                            l_ba = 0
                            l_q = 0
                            l_da = 0

                            #获取ref的碱基并涂色
                            color = 0
                            if line_31_arr[2] == 'a' or line_31_arr[2] == 'A':
                                color = 5 #灰色
                            elif line_31_arr[2] == 't' or line_31_arr[2] == 'T':
                                color = 6 
                            elif line_31_arr[2] == 'c' or line_31_arr[2] == 'C':
                                color = 7
                            else :
                                color = 8 


                            for i in range(b_len):
                                if l_b != b_len-1 and (bases[l_b] == ',' or bases[l_b] == '.' or bases[l_b] == 'n' or bases[l_b] == 'N' or bases[l_b] == '*'):
                                    l_b += 1
                                    if bases[l_b] == '+' or bases[l_b] == '-':
                                        bases_a[l_ba] = color
                                        bdfr_a[l_da] = 3
                                        l_da += 1
                                        l_ba += 1
                                        l_b += 1
                                        q = 0
                                        p = 1
                                        for j in range(l_b,b_len):
                                            if '0'<=bases[j]<='9':
                                                l_b += 1
                                            else: break
                                        for j in range(0,l_b):
                                            if '0'<=bases[l_b-1]<='9' :
                                                q = q + p*int(bases[l_b-1]) + 1
                                                p = p*10
                                                l_b -= 1
                                            else:   break
                                        l_b = l_b + q 
                                    elif bases[l_b-1] == '*': # *蓝色
                                        bases_a[l_ba] = color
                                        bdfr_a[l_da] = 3
                                        l_ba += 1
                                        l_da += 1
                                    else:
                                        bases_a[l_ba] = color
                                        bdfr_a[l_da] = 1
                                        l_ba += 1
                                        l_da += 1
                                elif l_b == b_len-1 and (bases[l_b] == ',' or bases[l_b] == '.' or bases[l_b] == 'n' or bases[l_b] == 'N'  or bases[l_b] == '*'):
                                    if bases[l_b] == '*':
                                        l_b += 1
                                        bases_a[l_ba] = color
                                        bdfr_a[l_da] = 3
                                        l_ba += 1
                                        l_da += 1
                                    else:
                                        l_b += 1
                                        bases_a[l_ba] = color
                                        bdfr_a[l_da] = 1
                                        l_ba += 1
                                        l_da += 1
                                elif bases[l_b] == '$':
                                    bases_a[l_ba] = 2
                                    bdfr_a[l_da] = 2
                                    l_da += 1
                                    l_ba += 1
                                    l_b += 1

                                elif bases[l_b] == 'a' or bases[l_b] == 'A':
                                    bases_a[l_ba] = 5 
                                    bdfr_a[l_da] = 4
                                    l_da += 1
                                    l_ba += 1
                                    l_b += 1
                                elif bases[l_b] == 't' or bases[l_b] == 'T':
                                    bases_a[l_ba] = 6 
                                    bdfr_a[l_da] = 4
                                    l_da += 1
                                    l_ba += 1
                                    l_b += 1
                                elif bases[l_b] == 'c' or bases[l_b] == 'C':
                                    bases_a[l_ba] = 7 
                                    bdfr_a[l_da] = 4
                                    l_da += 1
                                    l_ba += 1
                                    l_b += 1
                                elif bases[l_b] == 'g' or bases[l_b] == 'G':
                                    bases_a[l_ba] = 8 
                                    bdfr_a[l_da] = 4
                                    l_da += 1
                                    l_ba += 1
                                    l_b += 1
                                elif bases[l_b] == '^':
                                    l_b += 2

                                if l_q < q_len:
                                    if ord(qual[l_q]) >= 33 and ord(qual[l_q]) < 37:#1/2
                                        qual_a[l_q] = 9
                                        l_q += 1

                                    elif ord(qual[l_q]) >= 37 and ord(qual[l_q]) < 40:#1/4
                                        qual_a[l_q] = 8
                                        l_q += 1
                                    elif ord(qual[l_q]) >= 40 and ord(qual[l_q]) < 55:
                                        qual_a[l_q] = 6
                                        l_q += 1
                                    elif ord(qual[l_q]) >= 55 and ord(qual[l_q]) < 60:
                                        qual_a[l_q] = 5
                                        l_q += 1
                                    elif ord(qual[l_q]) >= 60 and ord(qual[l_q]) <= 126:
                                        qual_a[l_q] = 5
                                        l_q += 1

                                if l_b >= b_len:
                                    break

                            i = 0
                            for base in bases_a:
                                if base == 0 or (i == 32 and base != 2):
                                    break
                                elif i == 32 and base == 2:
                                    i -= 1
                                    if c == 30:
                                        for j in range(c,31):
                                            img_read_base[i][j] = 2
                                            img_base_quality[i][j] = 2
                                    else:
                                        img_read_base[i,c+1:] = 2
                                        img_base_quality[i,c+1:] = 2
                                        """ for j in range(c+1,31):
                                            img_read_base[i][j] = 2
                                            img_base_quality[i][j] = 2 """
                                    i += 1
                                else :    
                                    while img_read_base[i][c] == 2 and base != 2:
                                        i += 1
                                        if i == 32:
                                            break
                                    if i == 32:
                                        break
                                    if base == 2:
                                        if c == 30:
                                            for j in range(c,31):
                                                i_n = i - 1
                                                img_read_base[i_n][j] = 2
                                                img_base_quality[i_n][j] = 2
                                        else:
                                            i_n = i - 1
                                            img_read_base[i_n,c+1:] = 2
                                            img_base_quality[i_n,c+1:] = 2
                                            """ for j in range(c+1,31):
                                                i_n = i - 1
                                                img_read_base[i_n][j] = 2
                                                img_base_quality[i_n][j] = 2 """
                                    else:
                                        img_read_base[i][c] = base
                                        i += 1
                            #qual数组赋值
                            i = 0
                            for qu in qual_a:
                                if qu == 0 or (i == 32 and qu != 2):
                                    break
                                elif i == 32 and qu == 2:
                                    i -= 1
                                    if c == 30:
                                        for j in range(c,31):
                                            img_base_quality[i][j] = 2
                                    else:
                                        img_base_quality[i,c+1:] = 2
                                        """ for j in range(c+1,31):
                                            img_base_quality[i][j] = 2 """
                                    i += 1
                                else :    
                                    while img_base_quality[i][c] == 2 and qu != 2:
                                        i += 1
                                        if i == 32:
                                            break
                                    if i == 32:
                                        break
                                    if qu == 2:
                                        if c == 30:
                                            for j in range(c,31):
                                                i_n = i - 1
                                                img_base_quality[i_n][j] = 2
                                        else:
                                            i_n = i - 1
                                            img_base_quality[i_n,c+1:] = 2
                                            """ for j in range(c+1,31):
                                                i_n = i - 1
                                                img_base_quality[i_n][j] = 2 """
                                    else:
                                        img_base_quality[i][c] = qu
                                        i += 1    
                            #bdfr数组赋值
                            i = 0
                            for bdfr in bdfr_a:
                                if bdfr == 0 or (i == 32 and bdfr != 2):
                                    break
                                elif i == 32 and bdfr == 2:
                                    i -= 1
                                    if c == 30:
                                        for j in range(c,31):
                                            img_bdfr[i][j] = 2
                                    else:
                                        img_bdfr[i,c+1:] = 2
                                        """ for j in range(c+1,31):
                                            img_bdfr[i][j] = 2 """
                                    i += 1
                                else :    
                                    while img_bdfr[i][c] == 2 and bdfr != 2:
                                        i += 1
                                        if i == 32:
                                            break
                                    if i == 32:
                                        break
                                    if bdfr == 2:
                                        if c == 30:
                                            for j in range(c,31):
                                                i_n = i - 1
                                                img_bdfr[i_n][j] = 2
                                        else:
                                            i_n = i - 1
                                            img_bdfr[i_n,c+1:] = 2
                                            """ for j in range(c+1,31):
                                                i_n = i - 1
                                                img_bdfr[i_n][j] = 2 """
                                    else:
                                        img_bdfr[i][c] = bdfr
                                        i += 1
                            c += 1

                            if c == 31:

                                #对mapping数组、strand初始化，读取每一条reads的比对质量
                                #bamfile = pysam.AlignmentFile(path_bam, "rb")
                                k = 0
                                for bam in bam_array[jump:]: 
                                    bam_flag_2 = format(bam[0],'b')
                                    bam_len = len(bam_flag_2)
                                    bam_pass = 0
                                    if bam_len >= 9 :
                                        if bam_flag_2[bam_len-9] == '1':
                                            bam_pass = 1
                                    if bam_pass == 0:
                                        if k >= 32:
                                            break
                                        else :
                                            #初始化mapping
                                            if bam[1] >30:
                                                img_mq[k,:] = 5
                                                """ for j in range(0,31):
                                                    img_mq[k][j] = 5 """
                                            elif bam[1] >20:                                       
                                                img_mq[k,:] = 6
                                            elif bam[1] >10:
                                                img_mq[k,:] = 7
                                            elif bam[1] >=0:
                                                img_mq[k,:] = 8
                                            if bam_len >= 5:
                                                if bam_flag_2[bam_len-5] == '0' :
                                                    img_st[k,:] = 8
                                                    """ for j in range(0,31):
                                                        img_st[k][j] = 8 """
                                                elif bam_flag_2[bam_len-5] == '1':
                                                    img_st[k,:] = 5
                                                    """ for j in range(0,31):
                                                        img_st[k][j] = 5 """
                                            else:
                                                img_st[k,:] = 8
                                    k += 1
                                #bamfile.close()

                                for i in range(0,31):
                                    for j in range(0,32):
                                        if img_read_base[j][i] == 0 or img_read_base[j][i] == 2:
                                            img_mq[j][i] = img_st[j][i] = img_read_base[j][i]

                                img_rsv = array([[7 for x in range(31)] for y in range(32)],dtype=uint8)
                                #给img_rsv赋值
                                allele_count = 0
                                total_count = 0

                                for i in range(32):
                                    if img_bdfr[i][15] == 4:
                                        allele_count += 1
                                        total_count += 1
                                    elif  img_bdfr[i][15] != 0 and img_bdfr[i][15] != 2:
                                        total_count += 1

                                if total_count > 3 and allele_count/total_count > 1/4 or total_count <= 3 and allele_count > 1:

                                    for i in range(32):
                                        if img_bdfr[i][15] == 4:
                                            for j in range(31):
                                                if img_bdfr[i][j] == 0 or img_bdfr[i][j] == 2:
                                                    img_rsv[i][j] = img_bdfr[i][j]
                                                else :
                                                    img_rsv[i][j] = 5
                                        else :
                                            for j in range(31):
                                                if img_bdfr[i][j] == 0 or img_bdfr[i][j] == 2:
                                                    img_rsv[i][j] = img_bdfr[i][j]
                                else:
                                    for i in range(32):
                                        for j in range(31):
                                                if img_bdfr[i][j] == 0 or img_bdfr[i][j] == 2:
                                                    img_rsv[i][j] = img_bdfr[i][j]

                                """ #合并数组
                                for i in range(0,186):
                                    if i <= 30:
                                        for j in range(32):
                                            img_conbine[j][i] = img_read_base[j][i]
                                    elif i > 30 and i <= 61 :
                                        for j in range(32): 
                                            img_conbine[j][i] = img_base_quality[j][i-31]
                                    elif i > 61 and i <= 92 :
                                        for j in range(32): 
                                            img_conbine[j][i] = img_mq[j][i-62]
                                    elif i > 92 and i <= 123 :
                                        for j in range(32): 
                                            img_conbine[j][i] = img_st[j][i-93]
                                    elif i >  123 and i <= 154 :
                                        for j in range(32): 
                                            img_conbine[j][i] = img_rsv[j][i-124]
                                    else :
                                        for j in range(32): 
                                            img_conbine[j][i] = img_bdfr[j][i-155]

                                rd_img_array = []
                                for x in img_conbine.reshape(186*32):
                                    if x == 0:
                                        rd_img_array.append(255) # 白色-背景
                                    elif x == 1:
                                        rd_img_array.append(180) # 绿色-,.
                                    elif x == 2:
                                        rd_img_array.append(255) #白色-$
                                    elif x == 3:
                                        rd_img_array.append(40) # 蓝色-+ -
                                    elif x == 4:
                                        rd_img_array.append(20) # 红色-比错
                                    elif x == 5:
                                        rd_img_array.append(200) #灰色 边框
                                    elif x == 6:
                                        rd_img_array.append(149)
                                    elif x == 7:
                                        rd_img_array.append(100)
                                    elif x == 8:
                                        rd_img_array.append(50)
                                    elif x == 9:
                                        rd_img_array.append(20)
                                img = Image.new('L', (186, 32))
                                img.putdata(rd_img_array)
                                img.save(path_pic + '/' +ref + '-'  + 'snp' + flag_pos +  '.png') """
                                #输出数组

                                multi_channel_image = np.zeros((32, 31, 6), dtype=np.uint8)
                                multi_channel_image[:, :, 0] = img_read_base
                                multi_channel_image[:, :, 1] = img_base_quality
                                multi_channel_image[:, :, 2] = img_mq
                                multi_channel_image[:, :, 3] = img_st
                                multi_channel_image[:, :, 4] = img_rsv
                                multi_channel_image[:, :, 5] = img_bdfr

                                filename = path_pic + '/' + ref + '-' + 'snp' + flag_pos + '.npy'
                                np.save(filename,multi_channel_image)
                                #print(flag_pos,'已生成')

                                c = 0

                                #进行标注
                                true = True
                                false = False
                                null = None
                                """ if int(pos[tmp_pos]) == int(gt[tmp_gt]):#是snp
                                    num_snp += 1
                                    data={
                                    "version": "5.2.1",
                                    "flags": {
                                        "snp": true,
                                        "nosnp": false
                                    },
                                    "shapes": [],
                                    "imagePath": "..\\images\\snp" + str(int(pos[tmp_pos])) + ".png",
                                    "imageData": null,
                                    "imageHeight": 450,
                                    "imageWidth": 34
                                    }
                                    with open('../SNPTools/label/snp' + str(int(pos[tmp_pos])) + '.json', 'w') as f:
                                        json_str = json.dumps(data,indent=2)
                                        f.write(json_str)
                                    tmp_gt += 1

                                else :
                                    num_nosnp += 1
                                    data={
                                    "version": "5.2.1",
                                    "flags": {
                                        "snp": false,
                                        "nosnp": true
                                    },
                                    "shapes": [],
                                    "imagePath": "..\\images\\snp" + str(int(pos[tmp_pos])) + ".png",
                                    "imageData": null,
                                    "imageHeight": 450,
                                    "imageWidth": 34
                                    }
                                    with open('../SNPTools/label/snp' + str(int(pos[tmp_pos])) + '.json', 'w') as f:
                                        json_str = json.dumps(data,indent=2)
                                        f.write(json_str) """

                                """ while int(pos[tmp_pos]) > int(gt[tmp_gt]):
                                    tmp_gt += 1 """
                        tmp_pos += 1
                        #if count % 5000 == 0:
                            #print(flag_pos,'已生成')
                        #e = time()
                        #print('此点位用时',e-b)
                    """ while int(line_arr[1])+15 > int(pos[tmp_pos]):
                        tmp_pos += 1 """
                    if tmp_pos == len(pos) : #or tmp_gt == len(gt):   
                        break
                    """ while int(pos[tmp_pos]) > int(gt[tmp_gt]):
                        tmp_gt += 1
                        num_skipp += 1 """
                    #对其bam文件中比对质量和正反链    
                    bass = ''
                    for bas in basess:
                        if bas == '$' and bass != '^' :
                            jump += 1
                        bass = bas
                #print('snp数量',num_snp,'\n nosnp数量',num_nosnp)
                #break
    end_time = time()
    #print('总共用时',end_time - run_time)
    #print(num_skip,num_skipp)
