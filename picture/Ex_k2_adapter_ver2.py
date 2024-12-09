#阶梯式生成比对图片
#此版本将adapter变为四种颜色的碱基
#相较于ex-k2-adapter版本，此版本将base_a数组缩小
#相较于ver1，图片右边额外扩充以适配adapter被截断的情况
import pysam
from numpy import *
from PIL import Image
from time import *
from lxml import etree as ET
from itertools import islice

import os

#path_bam = 'chr20_d70.asm.bp.p_ctg.bam' #bam文件
#path_txt = 'chr20_d70.asm.bp.p_ctg.txt' #txt文件
#path_pic = 'pic/input/adapter' #图片存放路径
#path_label = 'adapter/label_0118' #标签存放地址

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


#out_file = open('base-ver1.txt','w') 

def adapter_pic(bam_dir, txt_add_dir, adapter_pic, chr):
    path_txt = txt_add_dir + '/' + chr + '_add.txt'
    path_bam = bam_dir + '/' + chr + '.bam'
    path_pic = adapter_pic
    #bamfile = pysam.AlignmentFile(bam, "rb")
    refs = []
    refs.append(chr)
    #bamfile.close()
    #print('Refrence Done!')
    ## 读取比对文件
    #begin_txt = time()
    """ txtfile = open(path_txt,"r") 
    lines = txtfile.readlines() """
    """ lines_len = 0
    with open(path_txt,'r') as lines:
        for line in lines:
           lines_len += 1 
    print('txt Done!')
    end_txt = time()
    run_txt = end_txt-begin_txt
    print ('txt读取时间：',run_txt) """
    a = 0
    #ref_n = 0
    bcf_n = 1
    num_ground = 0
    ## 以每一条参考read为一组
    #count = 0
    for ref in refs:
        count = 0
        lines_len = 0
        lines = list()
        begin_txt = time()
        with open(path_txt,'r') as f:
            for line in f:
                line_arr = line.split()
                if line_arr[0] == ref:
                    lines.append(line)
                    lines_len += 1
                elif lines_len != 0:
                    break
        print('txt Done!')
        end_txt = time()
        run_txt = end_txt-begin_txt
        print ('txt读取时间：',run_txt)
        #ref = 'chr20'
        ref_n = 0
        begin_r = time()
        print(ref)
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
        while flag == 0:
            ref_n += 1
            img_a = array([[0 for x in range(490)] for y in range(50)],dtype=uint8)      
            img_A = array([[0 for x in range(1471)] for y in range(901)],dtype=uint8)
            bases_a = array([0 for x in range(3000)])
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
            #with open(path_txt,'r') as lines:
            for line in lines[count_head:]:#islice(lines, count_head, lines_len):
                flag3 = 1 #当最后一条read的位置下面仍有read
                jump = 0#设置jump来修正txt文件中出现$而计数不对的问题
                if count_front >= 25:
                    flag = 0
                    print('next period',ref,int(line_arr[1]))
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
                    elif line_arr[0] != '':
                        ref = line_arr[0]
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
                                print('flas为',flag3)
                                pos = flag3
                                count_jump_head -= count_jump_head_tempory

                        ## 将比对信息写入数组
                        """ line1 = list(bases_a)
                        line1.insert(0,int(line_arr[1]))
                        line1.insert(0,',')
                        line1.insert(0,count_front)
                        line1.insert(0,',')
                        line1.insert(0,count_prior_10)
                        line1.insert(0,',')
                        line1.insert(0,count_jump)
                        line1.append('\n')
                        line = ''.join('%s' %id for id in line1)
                        out_file.writelines(line)  """
                        #tmp_indel = 0
                        i = 0
                        q_p = 0 #往前走
                        q_f = 0#往后走
                        basic = 0
                        for base in bases_a:
                            if (base == 0 and i >= 1) or i == 50:
                                break
                            elif base == 0:
                                pass
                            elif base < 0:
                                img_a[i][c] = 3
                                if img_a[i+1][c] == 3 or img_a[i+1][c-1] == 3:
                                    img_a[i][c] = 1
                                    i += 2
                                    pass
                                else :
                                    i += 1
                                    q_p = q_f = c
                                    basic = 0 - base
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
                                i += 2

                        c += 1



                if c == 445:
                    flag_save = 0
                    i_name = str(img_n)
                    r_name = str(ref_n)

                    filename='../images/' + 'ref' + r_name + '-' + i_name + '.png'
                    anno= GenAnnotations(filename, 1471, 901)#两类
                    #anno_adapter = GenAnnotations('../images-adapter/' + 'ref' + r_name + '-' + i_name + '.png', 1471, 901)
                    #anno_jiegou = GenAnnotations('../images-jiegou/' + 'ref' + r_name + '-' + i_name + '.png', 1471, 901)
                    for i in range(0,50):
                        qianhe = 0
                        adapter = 0
                        xmin = ymin = xmax = ymax = 0
                        for j in range(0,490):
                            if i % 2 == 0:#判断嵌合
                                """ if img_a[i][j] == 1 and qianhe == 0:
                                    qianhe = 1
                                    xmin = (j)*3
                                    ymin =(i)*3
                                if img_a[i][j] == 1 and qianhe != 0:
                                    qianhe += 1
                                if img_a[i][j] == 2 and qianhe < 1000:#符合条件
                                    xmax = (j)*3 + 3
                                    ymax =(i)*3 + 3
                                    anno.add_object("qianhe",xmin,ymin,xmax,ymax)
                                    anno.savefile('label/ref' + str(head) + '.xml')
                                    break
                                if img_a[i][j] == 2:
                                    break  """
                            else:#判断adapter

                                if img_a[i][j] == 3  and adapter == 0  :


                                    """ ada = 0
                                    for q in range(0,100):
                                        if j == 444:
                                            break
                                        if img_a[q][j] == 3  and q % 2 == 0:
                                            ada += 1
                                    if ada > 3 or ada == 0:
                                        adapter = 0
                                    else : """
                                    if j > 0:
                                        if img_a[i][j-1] != 3 and j == 1:

                                            adapter = 1

                                            ymin = (i-1)*18

                                            xmin = (j-1)*3
                                        elif img_a[i][j-1] != 3  and img_a[i][j-2] != 3  and j > 1:
                                            adapter = 1

                                            ymin = (i-1)*18

                                            xmin = (j-2)*3
                                        elif img_a[i][j-1] !=3  and img_a[i][j-2] == 3 and j > 1:
                                            adapter = 1

                                            ymin = (i-1)*18

                                            xmin = (j-1)*3
                                    elif j == 0:
                                            adapter = 1

                                            ymin = (i-1)*18

                                            xmin = (j)*3


                                elif img_a[i][j] == 3  and adapter != 0:
                                    adapter += 1


                                if (img_a[i][j] == 0 and 48 >= adapter >= 30 )or (j == 489 and 48 >= adapter >= 30) :#or img_a[i][j] == 3  and 48 >= adapter > 5:

                                    if j == 489 and adapter >= 30:
                                        xmax = (j)*3 + 3
                                    #elif img_a[i-1][j+1] == 3 and img_a[i][j] == 3 :
                                        #adapter = 0
                                    elif img_a[i][j] == 0 and adapter >= 30 :#or img_a[i][j] == 0 and j == adapter:#j == adapter 表示开头的adapter

                                        if j == 488:
                                            xmax = (j+1)*3 
                                        #elif j < 443 and img_a[i-1][j+1] == 3:
                                            #xmax = (j+1)*3
                                        else :
                                            xmax = (j+2)*3 


                                    ymax = i*18 + 18


                                    if xmax != 0  and ymax != 0 :

                                        flag_save = 1
                                        anno.add_object("adapter",xmin,ymin,xmax,ymax)
                                        num_ground += 1
                                        #anno.savefile(path_label + '/ref' + r_name + '-' + i_name + '.xml')
                                        adapter = 0
                                        #img.save(path_pic + '/' + 'ref' + r_name + '-' + i_name + '.png')
                                        xmin = ymin = xmax = ymax = 0
                                        """ if flag_adapter == 1 and adapter <=48: 
                                            anno.add_object("adapter",xmin,ymin,xmax,ymax)
                                            anno_adapter.add_object("adapter",xmin,ymin,xmax,ymax)
                                            num_ground += 1
                                            anno.savefile('data/label/ref' + r_name + '-' + i_name + '.xml')
                                            anno_adapter.savefile('data/label-adapter/ref' + r_name + '-' + i_name + '.xml')
                                            adapter = 0
                                            img.save(path_pic + '/' + 'ref' + r_name + '-' + i_name + '.png')
                                            img.save('data/images-adapter' + '/' + 'ref' + r_name + '-' + i_name + '.png')
                                            xmin = ymin = xmax = ymax = 0
                                        elif flag_jg == 1:
                                            anno.add_object("jiegou",xmin,ymin,xmax,ymax)
                                            anno_jiegou.add_object("jiegou",xmin,ymin,xmax,ymax)
                                            num_ground += 1
                                            anno.savefile('data/label/ref' + r_name + '-' + i_name + '.xml')
                                            anno_jiegou.savefile('data/label-jg/ref' + r_name + '-' + i_name + '.xml')
                                            adapter = 0
                                            img.save(path_pic + '/' + 'ref' + r_name + '-' + i_name + '.png')
                                            img.save('data/images-jg' + '/' + 'ref' + r_name + '-' + i_name + '.png')
                                            xmin = ymin = xmax = ymax = 0 """
                                if img_a[i][j] == 0 and adapter < 30 or img_a[i][j] == 0 and adapter > 48 :#or img_a[i][j] == 3 :
                                    if img_a[i][j] == 0 :
                                        adapter = 0
                                    """ elif img_a[i][j] == 3:
                                        if j == 444:
                                            pass 
                                        elif img_a[i-1][j+1] == 3 :
                                            adapter = 0 """
                    if flag_save > 0:
                        #当满足条件时才生成图片
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
                        img.save(path_pic + '/' + ref + '-' + r_name + '-' + i_name + '.png') 

                    c = 0               
                    img_n += 1
                    img_a = array([[0 for x in range(490)] for y in range(50)],dtype=uint8)
                    img_A = array([[0 for x in range(1471)] for y in range(901)],dtype=uint8)


                ## 每200个碱基生成一张图        
                #count_front为$计数，
                if count_front >= 25  or count == lines_len:
                    print(head,count_head,count_jump,count_jump_head)
                    flag_save = 0
                    i_name = str(img_n)
                    r_name = str(ref_n)

                    if count_front >= 25:
                        flag = 0


                    filename='../images/' + 'ref' + r_name + '-' + i_name + '-' + str(head) + '.png'
                    anno= GenAnnotations(filename, 1471, 901)
                    #anno_adapter = GenAnnotations('../images-adapter/' + 'ref' + r_name + '-' + i_name + '.png', 1471, 901)
                    #anno_jiegou = GenAnnotations('../images-jiegou/' + 'ref' + r_name + '-' + i_name + '.png', 1471, 901)
                    for i in range(0,50):
                        qianhe = 0
                        adapter = 0
                        xmin = ymin = xmax = ymax = 0
                        for j in range(0,490):
                            if i % 2 == 0:#判断嵌合
                                """ if img_a[i][j] == 1 and qianhe == 0:
                                    qianhe = 1
                                    xmin = (j)*3
                                    ymin =(i)*3
                                if img_a[i][j] == 1 and qianhe != 0:
                                    qianhe += 1
                                if img_a[i][j] == 2 and qianhe < 1000:#符合条件
                                    xmax = (j)*3 + 3
                                    ymax =(i)*3 + 3
                                    anno.add_object("qianhe",xmin,ymin,xmax,ymax)
                                    anno.savefile('label/ref' + str(head) + '.xml')
                                    break
                                if img_a[i][j] == 2:
                                    break  """
                            else:#判断adapter

                                if img_a[i][j] == 3  and adapter == 0  :



                                    """ ada = 0
                                    for q in range(0,100):
                                        if j == 444:
                                            break
                                        if img_a[q][j] == 3  and q % 2 == 0:
                                            ada += 1
                                    if ada > 3 or ada == 0:
                                        adapter = 0
                                    else : """
                                    if j > 0:
                                        if img_a[i][j-1] != 3 and j == 1:

                                            adapter = 1

                                            ymin = (i-1)*18

                                            xmin = (j-1)*3
                                        elif img_a[i][j-1] != 3  and img_a[i][j-2] != 3  and j > 1:
                                            adapter = 1

                                            ymin = (i-1)*18

                                            xmin = (j-2)*3
                                        elif img_a[i][j-1] !=3  and img_a[i][j-2] == 3 and j > 1:
                                            adapter = 1

                                            ymin = (i-1)*18

                                            xmin = (j-1)*3
                                    elif j == 0:
                                            adapter = 1

                                            ymin = (i-1)*18

                                            xmin = (j)*3


                                elif img_a[i][j] == 3  and adapter != 0:
                                    adapter += 1


                                if (img_a[i][j] == 0 and 48 >= adapter >= 30 )or (j == 489 and 48 >= adapter >= 30) :#or img_a[i][j] == 3  and 48 >= adapter > 5:

                                    if j == 489 and adapter >= 30:
                                        xmax = (j)*3 + 3
                                    #elif img_a[i-1][j+1] == 3 and img_a[i][j] == 3 :
                                        #adapter = 0
                                    elif img_a[i][j] == 0 and adapter >= 30 :#or img_a[i][j] == 0 and j == adapter:#j == adapter 表示开头的adapter

                                        if j == 488:
                                            xmax = (j+1)*3 
                                        #elif j < 443 and img_a[i-1][j+1] == 3:
                                            #xmax = (j+1)*3
                                        else :
                                            xmax = (j+2)*3 


                                    ymax =i*18 + 18


                                    if xmax != 0  and ymax != 0 :

                                        flag_save = 1
                                        anno.add_object("adapter",xmin,ymin,xmax,ymax)
                                        num_ground += 1
                                        #anno.savefile(path_label + '/ref' + r_name + '-' + i_name + '.xml')
                                        adapter = 0
                                        #img.save(path_pic + '/' + 'ref' + r_name + '-' + i_name + '.png')
                                        xmin = ymin = xmax = ymax = 0
                                        """ if flag_adapter == 1 and adapter <=48: 
                                            anno.add_object("adapter",xmin,ymin,xmax,ymax)
                                            anno_adapter.add_object("adapter",xmin,ymin,xmax,ymax)
                                            num_ground += 1
                                            anno.savefile('data/label/ref' + r_name + '-' + i_name + '.xml')
                                            anno_adapter.savefile('data/label-adapter/ref' + r_name + '-' + i_name + '.xml')
                                            adapter = 0
                                            img.save(path_pic + '/' + 'ref' + r_name + '-' + i_name + '.png')
                                            img.save('data/images-adapter' + '/' + 'ref' + r_name + '-' + i_name + '.png')
                                            xmin = ymin = xmax = ymax = 0
                                        elif flag_jg == 1:
                                            anno.add_object("jiegou",xmin,ymin,xmax,ymax)
                                            anno_jiegou.add_object("jiegou",xmin,ymin,xmax,ymax)
                                            num_ground += 1
                                            anno.savefile('data/label/ref' + r_name + '-' + i_name + '.xml')
                                            anno_jiegou.savefile('data/label-jg/ref' + r_name + '-' + i_name + '.xml')
                                            adapter = 0
                                            img.save(path_pic + '/' + 'ref' + r_name + '-' + i_name + '.png')
                                            img.save('data/images-jg' + '/' + 'ref' + r_name + '-' + i_name + '.png')
                                            xmin = ymin = xmax = ymax = 0 """
                                if img_a[i][j] == 0 and adapter < 30 or img_a[i][j] == 0 and adapter > 48 :#or img_a[i][j] == 3 :
                                    if img_a[i][j] == 0 :
                                        adapter = 0
                                    """ elif img_a[i][j] == 3:
                                        if j == 444:
                                            pass 
                                        elif img_a[i-1][j+1] == 3 :
                                            adapter = 0 """

                    if flag_save > 0:
                        #当满足条件时才生成图片
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
                        print ('图片生成时间：',run_time)
                        img.save(path_pic + '/' + ref + '-' + r_name + '-' + i_name + '.png')
                    img_n = 1
                    c = 0
            end_rr = time()
            run_rr = end_rr - begin_rr
            print('用时',run_rr)
        end_r = time()
        run_r = end_r-begin_r
        print ('本条read图片生成时间：',run_r)
        print('一共标注了adapter',num_ground)

