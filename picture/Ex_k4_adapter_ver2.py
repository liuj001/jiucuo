#阶梯式生成比对图片
#此版本将adapter变为四种颜色的碱基
#相较于ex-k2-adapter版本，此版本将base_a数组缩小
#相较于ver1，图片右边额外扩充以适配adapter被截断的情况
#相较于k2，在image_a数组赋值时便进行adapter标注，加快速度
#相较于K3，在read两边补充跳过的部分，以便识别adapter
import os
import pysam
import numpy as np
from numpy import *
from PIL import Image
from time import *
from lxml import etree as ET
from itertools import islice
class GenAnnotations:
    def __init__(self, filename, witdh, height, depth=3, foldername="VOC2007", sourcename="Unknown"):
        self.root = ET.Element("annotation")
        self.foleder = ET.SubElement(self.root, "folder")
        self.filename = ET.SubElement(self.root, "filename")
        self.source = ET.SubElement(self.root, "source")
        self.size = ET.SubElement(self.root, "size")
        self.width = ET.SubElement(self.size, "width")
        self.height = ET.SubElement(self.size, "height")
        self.depth = ET.SubElement(self.size, "depth")
        self.foleder.text = foldername
        self.filename.text = filename
        self.source.text = sourcename
        self.width.text = str(witdh)
        self.height.text = str(height)
        self.depth.text = str(depth)
    def savefile(self, filename):
        tree = ET.ElementTree(self.root)
        tree.write(filename, xml_declaration=False, encoding='utf-8', pretty_print=True)
    def add_object(self, label, xmin, ymin, xmax, ymax, tpose='Unspecified', ttruncated=0, tdifficult=0):
        object = ET.SubElement(self.root, "object")
        namen = ET.SubElement(object, "name")
        namen.text = label
        pose = ET.SubElement(object, "pose")
        pose.text = tpose
        truncated = ET.SubElement(object, "truncated")
        truncated.text = str(ttruncated)
        difficult = ET.SubElement(object, "difficult")
        difficult.text = str(tdifficult)
        bndbox = ET.SubElement(object, "bndbox")
        xminn = ET.SubElement(bndbox, "xmin")
        xminn.text = str(xmin)
        yminn = ET.SubElement(bndbox, "ymin")
        yminn.text = str(ymin)
        xmaxn = ET.SubElement(bndbox, "xmax")
        xmaxn.text = str(xmax)
        ymaxn = ET.SubElement(bndbox, "ymax")
        ymaxn.text = str(ymax)

#def images_make(ref):
def adapter_pic(bam_dir, txt_add_dir, adapter_dir, chr):
    #log_file = open(outdir+"/output.log", "a")
    #sys.stdout = log_file
    '''no_adapter = open('no_adapter.txt','r')
    no_adapter_list = no_adapter.readlines()
    no_adapter_list = [ ada.replace('\n','') for ada in no_adapter_list]'''
    ref = chr
    path_txt = txt_add_dir + '/' + ref + '_add.txt'
    path_pic = adapter_dir
    path_bam = bam_dir + '/' + ref + '.bam'
    #out_file = open('base-ver1.txt','w') 
    ## 获取reads的名称
    #print('Refrence Done!,adapter标注')
    ## 读取比对文件
    begin_txt = time()
    """ txtfile = open(path_txt,"r") 
    lines = txtfile.readlines() """
    end_txt = time()
    run_txt = end_txt-begin_txt
    #print ('txt读取时间：',run_txt)
    bcf_n = 1
    ##读取bam
    bamfile = pysam.AlignmentFile(path_bam, "rb")
    #out_txt1 = open('test_txt/' + ref + '-1.txt','w')
    #out_txt2 = open('test_txt/' + ref + '-2.txt','w')
    soft_clip_p = list()#前面跳过部分
    soft_clip_f = list()
    tmp_b = 0
    for bam in bamfile:
        # if bam.reference_name != ref and tmp_b > 0:
        #     break
        # if bam.reference_name == ref:
            # if tmp_b <= 22575 and ref == 'utg0002050l':
            #     pass
        
        if tmp_b >= 0:
            if bam.cigar[0][0] == 4:
                '''if (bam.query_name in no_adapter_list):
                    no_adapter_list.remove(bam.query_name'''
                if bam.cigar[0][1] < 30:
                    pass
                elif 45 >= bam.cigar[0][1] >= 30:
                    soft_clip_p.append([tmp_b,bam.cigar[0][1]])
                else:
                    soft_clip_p.append([tmp_b,45])
            if bam.cigar[-1][0] ==  4 :
                '''if (bam.query_name in no_adapter_list):
                    no_adapter_list.remove(bam.query_name)'''
                if bam.cigar[-1][1] < 30:
                    pass
                elif 45 >= bam.cigar[-1][1] >= 30:
                    soft_clip_f.append([tmp_b,bam.cigar[-1][1]])
                else:
                    soft_clip_f.append([tmp_b,45])
        tmp_b += 1
    #print(len(soft_clip_p))
    #print(len(soft_clip_f))
    #for row in soft_clip_p:
        #out_txt1.write(' '.join(map(str, row)) + '\n')
    #out_txt1.close()
    #for row in soft_clip_f:
        #out_txt2.write(' '.join(map(str, row)) + '\n')
    #out_txt2.close()
    bamfile.close()
    #print('bam done')
    ## 以每一条参考read为一组
    '''in_file = open(path_txt,'r')
    lines = in_file.readlines()'''
    lines_len = 0
    with open(path_txt,'r') as lines:
        for line in lines:
           lines_len += 1 
    #print(ref + 'txt Done!')
    num_ground = 0
    count = 0
    ref_n = 0
    begin_r = time()
    #print(ref)
    count_prior = 0
    count_prior_10 = 0
    head = 0
    img_a = array([[0 for x in range(490)] for y in range(50)],dtype=uint8)
    img_A = array([[0 for x in range(1471)] for y in range(901)],dtype=uint8) 
    flag = 0
    flag2 = 1
    c = 0
    pos = 1
    pos_f = 1
    count_jump_head = 0 #记录每次head之前的count_jump数，防止重复遍历
    count_head = 0 #记录head所在的行数
    flag_save = 0
    while flag == 0:
        if len(soft_clip_p) == 0 and len(soft_clip_f) == 0:
            #print('have no prior or front adapter ')
            # break
            pass
        ref_n += 1
        img_a = array([[0 for x in range(490)] for y in range(50)],dtype=uint8)      
        img_A = array([[0 for x in range(1471)] for y in range(901)],dtype=uint8)
        bases_a = array([0 for x in range(50000)])
        flag = 1
        flag1 = 0
        flag2 = 0
        count_front = 0
        count_jump = count_jump_head #count_jump代表每次循环之前结束的read数量
        count = count_head
        count_prior_10 = count_prior
        img_n = 1
        begin_rr = time()
        if pos > 1:
            head -= 1
            count_head -= 1 
        with open(path_txt,'r') as lines:
            for line in islice(lines, count_head, lines_len):
                flag3 = 1 #当最后一条read的位置下面仍有read
                jump = 0#设置jump来修正txt文件中出现$而计数不对的问题
                if count_front >= 25:
                    flag = 0
                    #print('next period',int(line_arr[1]),ref)
                    break
                if flag1 == 0:
                    count_head += 1
                count += 1
                line_arr = line.split()
    
                if int(line_arr[1]) <= head :
                    pass
                    """ p_p = ''
                    for p in line:
                        if p == '$' and p_p != '^':
                            count_jump += 1
                        p_p = p  """
                            
                else :
                    num = 0 #num的作用：如果当前开始（^）的位置有$，不计入count_jump，后续遍历中会计入count_jump
                    p_p = ''
                    for p in line_arr[4]:
                        if flag2 == 1:
                            break
                        if p == '$'and p_p != '^':
                            count_jump_head += 1
                            count_jump += 1
                            num += 1
                        if p == '^' and flag2 == 0:
                            flag2 = 1
                            count_jump -= num
                            count_jump_head -= num 
                            if pos > 1:
                                pos_f = pos
                                pos = 1
                        p_p = p
    
                    if flag2 == 0:
                        pass
                # print(line_arr)
                ## 比对到当前参考read的reads,提取比对信息
                    elif line_arr[0] == ref or 1 == 1:
                        
                        tmp_indel = 0
                        bases = list(line_arr[4])
                        b_len = len(bases)
                        #indel = array([0 for x in range(b_len)],dtype=uint8) #记录插入碱基的四种不同颜色
                        # print(b_len)
                        l_b = 0
                        l_ba = 0
                        count_jump_temporary = 0
                        count_jump_head_tempory = 0
                        ## 将比对信息转换成相对应的数字
                        #检测到,.后无+-，bases_a为1，有+-为3
                        #检测到$，bases_a为2
                        #检测到atcg为4
                        #其余全为0
                        for i in range(b_len):
                            ## 控制bases_a的长度 增加速度
                            if (l_ba + count_jump - count_prior_10 > 60 or l_b < 0) and bases[l_b] != '^' and bases[l_b-1] != '^' and bases[l_b-2] != '^':
                                break
                            
                            #count_jump为跳过计数
                            if l_b != b_len-1 and (bases[l_b] == ',' or bases[l_b] == '.' or bases[l_b] == '*' or bases[l_b] == 'a' or bases[l_b] == 'A'or bases[l_b] == 'c' or bases[l_b] == 'C' or bases[l_b] == 'g' or bases[l_b] == 'G' or bases[l_b] == 't' or bases[l_b] == 'T' or bases[l_b] == 'N' or bases[l_b] == 'n'):
                                l_b += 1
                                if bases[l_b] == '+' or bases[l_b] == '-':
                                    while bases_a[l_ba + count_jump - count_prior_10 ] == 2:
                                        l_ba += 1
                                    if l_ba >= count_prior_10 - count_jump and bases[l_b] == '-':
                                        bases_a[l_ba + count_jump - count_prior_10 ] = 1
                                    elif l_ba >= count_prior_10 - count_jump and bases[l_b] == '+':
                                        bases_a[l_ba + count_jump - count_prior_10] = -1
                                    l_b += 1
                                    l_ba += 1
                                    q = 0
                                    p = 1
                                    cc = 0
                                    for j in range(l_b,b_len):
                                        if '0'<=bases[j]<='9':
                                            l_b += 1
                                        else: break
                                    l_bb = l_b #用l_b记住此位置，然后向后遍历插入碱基
                                    for j in range(0,l_b):
                                        if '0'<=bases[l_b-1]<='9' :
                                            q = q + p*int(bases[l_b-1]) + 1
                                            p = p*10
                                            l_b -= 1
                                            cc += 1
                                        else:   break
                                    l_b = l_b + q 
    
                                    if bases_a[l_ba - 1 + count_jump - count_prior_10] == -1:
                                        if q - cc > 46:
                                            bases_a[l_ba - 1 + count_jump - count_prior_10] = -46
                                        else :
                                            bases_a[l_ba - 1 + count_jump - count_prior_10] = -q + cc #此位置置位蓝色，后续位置移到下一行，并记录碱基颜色
                                        """ for j in range(0,q-cc):
                                            if bases[l_bb] == 'a' or bases[l_bb] == 'A':
                                                indel[tmp_indel] = 6
                                            elif bases[l_bb] == 't' or bases[l_bb] == 'T':
                                                indel[tmp_indel] = 7
                                            elif bases[l_bb] == 'c' or bases[l_bb] == 'C':
                                                indel[tmp_indel] = 8
                                            elif bases[l_bb] == 'g' or bases[l_bb] == 'G':
                                                indel[tmp_indel] = 9
                                            tmp_indel += 1
                                            l_bb += 1 """
                                        
                                elif bases[l_b-1] == '*': # *为蓝色
                                    while bases_a[l_ba+ count_jump - count_prior_10] == 2:
                                        l_ba += 1
                                    if l_ba >= count_prior_10 - count_jump:
                                        bases_a[l_ba +count_jump - count_prior_10] = 7
                                    l_ba += 1
                                elif bases[l_b-1] == ',' or bases[l_b-1] == '.':
                                    while bases_a[l_ba + count_jump - count_prior_10] == 2:
                                        l_ba += 1                           
                                    if l_ba >= count_prior_10 -count_jump:
                                        bases_a[l_ba + count_jump - count_prior_10] = 1
                                    l_ba += 1
                                else :
                                    while bases_a[l_ba + count_jump- count_prior_10] == 2:
                                        l_ba += 1                           
                                    if l_ba >= count_prior_10 -count_jump:
                                        bases_a[l_ba + count_jump - count_prior_10 ] = 4
                                    l_ba += 1
                            elif l_b == b_len-1 and (bases[l_b] == ',' or bases[l_b] == '.' or bases[l_b] == '*' or bases[l_b] == 'a' or bases[l_b] == 'A'or bases[l_b] == 'c' or bases[l_b] == 'C' or bases[l_b] == 'g' or bases[l_b] == 'G' or bases[l_b] == 't' or bases[l_b] == 'T' or bases[l_b] == 'N' or bases[l_b] == 'n'):
                                if bases[l_b] == '*':
                                    l_b += 1
                                    while bases_a[l_ba + count_jump - count_prior_10] == 2:
                                        l_ba += 1                            
                                    if l_ba >= count_prior_10 - count_jump:
                                        bases_a[l_ba +count_jump - count_prior_10] = 7
                                    l_ba += 1
                                elif bases[l_b] == ',' or bases[l_b] == '.':
                                    l_b += 1
                                    while bases_a[l_ba + count_jump - count_prior_10] == 2:
                                        l_ba += 1                            
                                    if l_ba >= count_prior_10 -count_jump:
                                        bases_a[l_ba +count_jump - count_prior_10] = 1
                                    l_ba += 1
                                else :
                                    l_b += 1
                                    while bases_a[l_ba + count_jump - count_prior_10] == 2:
                                        l_ba += 1                            
                                    if l_ba >= count_prior_10 -count_jump:
                                        bases_a[l_ba +count_jump - count_prior_10] = 4
                                    l_ba += 1
                            elif bases[l_b] == '$' and bases[l_b-1] != '^' :#如果bases[l_b] == '$'，则上一行肯定有值，bases_a[l_ba] == 2 跳行需要先跳到下一行有值的行，在剪去跳数
                                if bases_a[l_ba + count_jump - count_prior_10] == 2:
                                    l_ba1 = 1
                                    while bases_a[l_ba + count_jump - count_prior_10] == 2:
                                        l_ba += 1
                                        l_ba1 += 1
                                    if count_prior_10 + 24 - count_jump>= l_ba - l_ba1 >= count_prior_10 - count_jump:
                                        count_front += 1
                                    if l_ba >= count_prior_10 +1 - count_jump:
                                        bases_a[l_ba - l_ba1 + count_jump - count_prior_10] = 2
                                    else :
                                        count_jump_temporary += 1 #此时count_jump不要+，
                                        jump += 1
                                    if flag1 == 0:#说明未到达head count_jump_head_tempory+1,此时不管结束的位置是否在图的区间内，因为count_jump_head_tempory是全局jump，所有都要加1
                                            count_jump_head_tempory += 1
                                else :
                                    if count_prior_10 + 24 - count_jump >= l_ba - 1 >= count_prior_10 - count_jump:
                                        count_front += 1          
                                    if l_ba >= count_prior_10 + 1 - count_jump:
                                        bases_a[l_ba -1 + count_jump- count_prior_10 ] = 2
                                    else :
                                        count_jump_temporary += 1
                                        jump += 1
                                    if flag1 == 0:#说明未到达head
                                        count_jump_head_tempory += 1
                                l_b += 1
                            #elif bases[l_b] == 'a' or bases[l_b] == 'A'or bases[l_b] == 'c' or bases[l_b] == 'C' or bases[l_b] == 'g' or bases[l_b] == 'G' or bases[l_b] == 't' or bases[l_b] == 'T' or bases[l_b] == 'N' or bases[l_b] == 'n':
                            #    while bases_a[l_ba + count_jump] == 2:
                            #        l_ba += 1
                            #    if l_ba >= count_prior_10 -count_jump:
                            #        bases_a[l_ba +count_jump] = 4 
                            #    l_b += 1
                            #    l_ba += 1 
                            elif bases[l_b] == '^':  
                                if pos_f > 1:#pos会导致count_prior计数不对 pos_f>1证明起点比对
                                    l_b = b_len - 3*(pos_f - 1)
                                    #调整l_ba
                                    if l_ba == 0:
                                        l_ba += count_prior_10
                                    #l_ba += pos_f - 1
                                    if  flag1 == 0 and count_prior - count_prior_10 < 25 :
                                        count_prior += 1 
                                    elif count_prior - count_prior_10 == 25 :
                                        flag3 += 1 #该位置增加了多条read
                                    pos_f -= 1
                                    l_b += 2
    
                                else:
                                    l_b += 2
                                    if  flag1 == 0 and count_prior - count_prior_10 < 25 :
                                        count_prior += 1
                                    elif count_prior - count_prior_10 == 25 :
                                        flag3 += 1 #该位置增加了多条read
                            if l_b >= b_len:
                                break
    
                        count_jump += count_jump_temporary
                        count_jump_head += count_jump_head_tempory
                        # print(bases)
                        # print(bases_a)
                        ''' if jump >= 1 and count_prior %25 != 0:#说明出现了$
                            q = 0
                            while bases_a[q] == 0:
                                q += 1
                            while bases_a[q] != 0:
                                q += 1
                            while jump != 0:
                                bases_a[q-1] =  0
                                q -= 1
                                jump -= 1 '''
                        
                        #count_prior为^计数，当增到10的时候令flag1为1，便不再增长
                        if count_prior % 25 == 0 and count_prior != 0 and flag1 == 0:
                            head = int(line_arr[1])
                            flag1 = 1
                            if flag3 > 1:#说明结束的位置仍有多条新加的read
                                #print('flas为',flag3)
                                pos = flag3
                                count_jump_head -= count_jump_head_tempory
                            
                        ## 将比对信息写入数组
                        
                        '''line1 = list(bases_a)
                        line1.insert(0,int(line_arr[1]))
                        line1.insert(0,',')
                        line1.insert(0,count_front)
                        line1.insert(0,',')
                        line1.insert(0,count_prior_10)
                        line1.insert(0,',')
                        line1.insert(0,count_jump)
                        line1.append('\n')
                        line = ''.join('%s' %id for id in line1)
                        out_file.writelines(line) ''' 
                        #tmp_indel = 0
                        if c == 0:
                            filename='../images/' + ref + '-' + str(ref_n) + '-' + str(img_n) + '.png'
                            anno= GenAnnotations(filename, 1471, 901)#两类
                        i = 0
                        q_p = 0 #往前走
                        q_f = 0#往后走
                        basic = 0
                        for base in bases_a:                       
                            if (base == 0 and i >= 1) or i == 50:
                                break
                            elif base == 0:
                                pass
                            elif base < 0: #插入
                                img_a[i][c] = 3
                                if img_a[i+1][c] == 3 or img_a[i+1][c-1] == 3:
                                    img_a[i][c] = 1
                                    i += 2
                                    pass
                                else :
                                    basic = 0 - base
                                    if 48 >= basic >= 30:#进行标注
                                        num_ground += 1
                                        flag_save = 1
                                        '''xmin = ymin = xmax = ymax = 0
                                        if c > 0:
                                            if img_a[i][c-1] != 3 and c == 1: 
                                                ymin = (i)*18                                           
                                                xmin = (c-1)*3
                                            elif img_a[i][c-1] != 3  and img_a[i][c-2] != 3  and c > 1:                                          
                                                ymin = (i)*18                                          
                                                xmin = (c-2)*3
                                            elif img_a[i][c-1] !=3  and img_a[i][c-2] == 3 and c > 1: 
                                                ymin = (i)*18                                         
                                                xmin = (c-1)*3
                                        elif c == 0:
                                                ymin = (i)*18
                                                xmin = (c)*3
                                        if c + basic >= 489:
                                            xmax = (489)*3 + 3
                                        elif c + basic == 488:
                                            xmax = (c + basic + 1)*3
                                        else:
                                            xmax = (c + basic + 2)*3
                                        ymax = (i+1)*18 + 18
                                        anno.add_object("adapter",xmin,ymin,xmax,ymax)'''

                                    i += 1
                                    q_p = q_f = c
                                    for bas in range(0,basic):
                                        if q_f != 490:
                                            img_a[i][q_f] = 3
                                            #tmp_indel += 1
                                            q_f += 1
                                        else :
                                            break
                                            """ img_a[i][q_p-1] = indel[tmp_indel]
                                            q_p -= 1
                                            tmp_indel += 1
                                            if q_p == 0:
                                                break """
                                    i += 1
                            else :  
                                img_a[i][c] = base
                                ##k4增加部分，在此加入左右跳过部分
                                #左边加入
                                if len(soft_clip_p) != 0:
                                    if 24 + count_prior_10 < soft_clip_p[0][0]:
                                        pass
                                    elif (i/2 + count_prior_10 == soft_clip_p[0][0]) and (c == 0 or img_a[i][c-1] == 0) and base == 1:
                                        img_a[i][c] = 3
                                        i += 1
                                        img_a[i][c:c+soft_clip_p[0][1]] = 3
                                        i -= 1
                                        del soft_clip_p[0]
                                        flag_save = 1
                                #右边加入
                                if len(soft_clip_f) != 0:
                                    for f,clip in enumerate(soft_clip_f):
                                        if clip[0] > count_prior_10 + 24:   
                                            break
                                        elif i/2 + count_prior_10 == clip[0] and img_a[i][c-1] == 1 and base == 2:
                                            
                                            img_a[i][c-1] = 3
                                            i += 1
                                            img_a[i][c-1:c + clip[1] - 1] = 3
                                            i -= 1
                                            del soft_clip_f[f]
                                            flag_save = 1
                                        
                                ##   
                                i += 2
                                                            
                        c += 1
                        
                        
                
                if c == 445:
                    i_name = str(img_n)
                    r_name = str(ref_n)
                    
                    if flag_save > 0:
                        flag_save = 0
                        #anno.savefile('label_adapter/'+ ref + '-' + r_name + '-' + i_name + '.xml')
                        #当满足条件时才生成图片
                        ##保持numpy
                        i1 = 0
                        j1 = 0
                        #为445*100的图补边框->1999*501
                        for i in range(0,50): 
                            for j in range(0,490):
                                if img_a[i][j] != 0 and img_a[i][j] != 2:
                                    j1 = 3*j
                                    i1 = 18*i
                                    img_A[i1:i1+18,j1] = 5
                                    """ for k in range(0,18):
                                            img_A[i1 + k][j1] = 5 """
                                    img_A[i1,j1:j1 + 3] = 5
                                    """ for k in range(0,3):
                                            img_A[i1][j1 + k] = 5 """
                                    img_A[i1+1:i1 + 18,j1+1:j1 + 3] = img_a[i][j]
                                    """ for k in range(1,18):
                                            for g in range(1,3):
                                                img_A[i1 + k][j1 + g] = img_a[i][j] """
                                    #img_A[i1][j1] = img_A[i1+1][j1] = img_A[i1+2][j1] = img_A[i1][j1+1] = img_A[i1][j1+2] = 5
                                    #img_A[i1+1][j1+1] = img_A[i1+1][j1+2] = img_A[i1+2][j1+1] = img_A[i1+2][j1+2] = img_a[i][j]
                                    if j == 489:
                                        img_A[18*i:18*i+18,1470]  = 5
                                        """ for k in range(0,18):
                                            img_A[18*i+k][1470]  = 5 """
                                        #img_A[4*i][1470] = img_A[4*i+1][1470] = img_A[4*i+2][1470]  = 5
                                else :
                                    j1 = 3*j
                                    i1 = 18*i
                                    """ for k in range(0,18):
                                        for g in range(0,3):
                                            img_A[i1 + k][j1 + g] = img_a[i][j] """
                                    #img_A[i1][j1] = img_A[i1+1][j1] = img_A[i1+2][j1] = img_A[i1][j1+1] = img_A[i1][j1+2] = img_a[i][j] 
                                    #img_A[i1+1][j1+1] = img_A[i1+1][j1+2] = img_A[i1+2][j1+1] = img_A[i1+2][j1+2] = img_a[i][j]
                                    if img_a[i][j-1] != 0 and img_a[i][j-1] != 2 and j>=1:
                                        j1 = 3*j
                                        i1 = 18*i
                                        img_A[i1:i1 + 18,j1] = 5
                                        """ for k in range(0,18):
                                            img_A[i1 + k][j1] = 5 """
                                        #img_A[i1][j1] = img_A[i1+1][j1] = img_A[i1+2][j1] = 5
                                    if img_a[i-1][j] != 0 and img_a[i-1][j] != 2 and i >= 1:
                                        j1 = 3*j
                                        i1 = 18*i
                                        img_A[i1,j1:j1 + 3] = 5
                                        """ for k in range(0,3):
                                            img_A[i1][j1 + k] = 5 """
                                        #img_A[i1][j1] = img_A[i1][j1+1] = img_A[i1][j1+2] = 5
                            if i == 49 :
                                for j in range(0,490):
                                    if img_a[i][j] != 2 and img_a[i][j] != 0:
                                        img_A[18*i+18,3*j:3*j + 3] = 5
                                        """ for k in range(0,3):
                                            img_A[18*i+18][3*j + k] = 5 """
                                        #img_A[6*i+6][6*j] = img_A[6*i+6][6*j+1] = img_A[6*i+6][3*j+2] = 5
                            
                        # print(img_a)
                        begin_time = time()
    
                        ## 转换为彩色
                        img_array = []
                        for x in img_A.reshape(1471*901):
                            if x == 0:
                                img_array.append((255,255,255)) # 白色-背景
                            elif x == 1:
                                img_array.append((0,255,0)) # 绿色-,.
                            elif x == 2:
                                img_array.append((255,255,255)) #白色-$
                            elif x == 3:
                                img_array.append((0,0,255)) # 蓝色-+ -
                            elif x == 4:
                                img_array.append((255,0,0)) # 红色-比错
                            elif x == 5:
                                img_array.append((0,0,0)) #黑色 边框
                            elif x == 7:
                                img_array.append((0,255,255))
    
                        ## 生成图片
                        img = Image.new('RGB', (1471, 901))
                        img.putdata(img_array)
                        #img.save(path_pic + '/' + 'ref' + r_name + '-' + i_name + '.png') 
                        img.save(path_pic + '/' + ref + '-' + r_name + '-' + i_name + '.png')
                        
                    c = 0               
                    img_n += 1
                    img_a = array([[0 for x in range(490)] for y in range(50)],dtype=uint8)
                    img_A = array([[0 for x in range(1471)] for y in range(901)],dtype=uint8)
    
    
                ## 每200个碱基生成一张图        
                #count_front为$计数，
                if count_front >= 25  or count == lines_len:
                    #print(head,count_head,count_jump,count_jump_head)
                    i_name = str(img_n)
                    r_name = str(ref_n)
                    
                    if count_front >= 25:
                        flag = 0
    

                    if flag_save > 0:
                        flag_save = 0
                        #anno.savefile('label_adapter/'+ ref + '-' + r_name + '-' + i_name + '.xml')
                        #当满足条件时才生成图片
                        ##保存numpy
                        i1 = 0
                        j1 = 0
                        #为445*100的图补边框->1999*501
                        for i in range(0,50): 
                            for j in range(0,490):
                                if img_a[i][j] != 0 and img_a[i][j] != 2:
                                    j1 = 3*j
                                    i1 = 18*i
                                    img_A[i1:i1+18,j1] = 5
                                    """ for k in range(0,18):
                                            img_A[i1 + k][j1] = 5 """
                                    img_A[i1,j1:j1 + 3] = 5
                                    """ for k in range(0,3):
                                            img_A[i1][j1 + k] = 5 """
                                    img_A[i1+1:i1 + 18,j1+1:j1 + 3] = img_a[i][j]
                                    """ for k in range(1,18):
                                            for g in range(1,3):
                                                img_A[i1 + k][j1 + g] = img_a[i][j] """
                                    #img_A[i1][j1] = img_A[i1+1][j1] = img_A[i1+2][j1] = img_A[i1][j1+1] = img_A[i1][j1+2] = 5
                                    #img_A[i1+1][j1+1] = img_A[i1+1][j1+2] = img_A[i1+2][j1+1] = img_A[i1+2][j1+2] = img_a[i][j]
                                    if j == 489:
                                        img_A[18*i:18*i+18,1470]  = 5
                                        """ for k in range(0,18):
                                            img_A[18*i+k][1470]  = 5 """
                                        #img_A[4*i][1470] = img_A[4*i+1][1470] = img_A[4*i+2][1470]  = 5
                                else :
                                    j1 = 3*j
                                    i1 = 18*i
                                    """ for k in range(0,18):
                                        for g in range(0,3):
                                            img_A[i1 + k][j1 + g] = img_a[i][j] """
                                    #img_A[i1][j1] = img_A[i1+1][j1] = img_A[i1+2][j1] = img_A[i1][j1+1] = img_A[i1][j1+2] = img_a[i][j] 
                                    #img_A[i1+1][j1+1] = img_A[i1+1][j1+2] = img_A[i1+2][j1+1] = img_A[i1+2][j1+2] = img_a[i][j]
                                    if img_a[i][j-1] != 0 and img_a[i][j-1] != 2 and j>=1:
                                        j1 = 3*j
                                        i1 = 18*i
                                        img_A[i1:i1 + 18,j1] = 5
                                        """ for k in range(0,18):
                                            img_A[i1 + k][j1] = 5 """
                                        #img_A[i1][j1] = img_A[i1+1][j1] = img_A[i1+2][j1] = 5
                                    if img_a[i-1][j] != 0 and img_a[i-1][j] != 2 and i >= 1:
                                        j1 = 3*j
                                        i1 = 18*i
                                        img_A[i1,j1:j1 + 3] = 5
                                        """ for k in range(0,3):
                                            img_A[i1][j1 + k] = 5 """
                                        #img_A[i1][j1] = img_A[i1][j1+1] = img_A[i1][j1+2] = 5
                            if i == 49 :
                                for j in range(0,490):
                                    if img_a[i][j] != 2 and img_a[i][j] != 0:
                                        img_A[18*i+18,3*j:3*j + 3] = 5
                                        """ for k in range(0,3):
                                            img_A[18*i+18][3*j + k] = 5 """
                                        #img_A[6*i+6][6*j] = img_A[6*i+6][6*j+1] = img_A[6*i+6][3*j+2] = 5
                            
                        # print(img_a)
                        begin_time = time()
    
                        ## 转换为彩色
                        img_array = []
                        for x in img_A.reshape(1471*901):
                            if x == 0:
                                img_array.append((255,255,255)) # 白色-背景
                            elif x == 1:
                                img_array.append((0,255,0)) # 绿色-,.
                            elif x == 2:
                                img_array.append((255,255,255)) #白色-$
                            elif x == 3:
                                img_array.append((0,0,255)) # 蓝色-+ -
                            elif x == 4:
                                img_array.append((255,0,0)) # 红色-比错
                            elif x == 5:
                                img_array.append((0,0,0)) #黑色 边框
                            elif x == 7:
                                img_array.append((0,255,255))
    
    
                        ## 生成图片
                        img = Image.new('RGB', (1471, 901))
                        img.putdata(img_array)
                        end_time = time()
                        run_time = end_time-begin_time
                        #print ('图片生成时间：',run_time)
                        #img.save(path_pic + '/' + 'ref' + r_name + '-' + i_name + '.png')
                        img.save(path_pic + '/' + ref + '-' + r_name + '-' + i_name + '.png')
                        
                    img_n = 1
                    c = 0
            end_rr = time()
            run_rr = end_rr - begin_rr
            #print('用时',run_rr)
    end_r = time()
    run_r = end_r-begin_r
    #print ('本条read图片生成时间：',run_r)
    #print('一共标注了adapter',num_ground)
    #out_txt_done = open('test_txt/' + ref + '-done.txt','w')
    #out_txt_done.close()
    #in_file.close()
    #log_file.close()
    #sys.stdout = sys.__stdout__ 