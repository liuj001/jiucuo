import os
import pysam
from itertools import islice
import time
import gzip

#此版本固定了删除和插入碱基个数为5
#此版本使用多个txt文件修改错误 v3
#此版本 将找read名字的操作删除 lines_read_name lines_read_name_tmp v4
#此版本 将read和bamlines按照txt中的深度和点位存入内存 而非一次性 v5
def correct(bam_dir, txt_dir, ec_dir, c_reads_dir, chr):
    
    #directory = '/root/autodl-tmp/test/' #txt文件和ec文件目录
    path_bam = bam_dir + '/' + chr + '.bam'  #bam文件
    #in_read = gzip.open('Z.mays/CRR302668.fastq.gz','rt') #初始read
    out_read = gzip.open(c_reads_dir+'/'+chr+'.ec.fastq.gz','wt')  #修改后的read
    #out_check = open('check.txt','w')
    ec_list = list()
    txt_list = list()
    ecf = chr+'_ec.txt'
    txtf = chr+'.txt'
    ec_list.append(ecf)
    txt_list.append(txtf)
    # for file in sorted(os.listdir(directory)):#获得txt文件和ec文件列表
    #     if os.path.isfile(os.path.join(directory, file)) and os.path.splitext(file)[1] == '.txt':#获取后缀为txt的文件名
    #         if 'ec' in file:#获取ec文件
    #             ec_list.append(file)
    #         elif 'bcf' not in file:#获取txt文件
    #             txt_list.append(file)
    #print(txt_list,ec_list)
    #if len(txt_list) != len(ec_list):
        #print('文件数量不对应！')
    time_fs = time.time()
    #lines_txt = txt_file.readlines()
    #print('Correction Begin')
    #prior_temporary = 0
    #prior_num = 0
    #ci_temporary = 0
    time_100000_b = 0
    time_100000_e = 0
    #print('--------------------------len(txt_list):',len(txt_list))

    for k in range(len(txt_list)):
        
        #print(ec_list[k])
        #print(txt_list[k])

        pos = list()
        pos_tail = 0 #防止最后一个元素在Pos中弹出后在进入是0
        num_pos = 1

        #print('No.',k+1,'ptg===========')
        lines_read = list()
        bamlines = list()#储存bam信息，防止循环遍历
        #print('-------------------------bam:',path_bam)
        bamfile = pysam.AlignmentFile(path_bam,'rb')
        refs = list()
        refs.append(chr)
        #print('-------------------------refs:',refs)

        len_bam = 0
        for bam in bamfile:
            if bam.reference_name == refs[k]:
                len_bam += 1
                """ bamlines.append(bam)
                lines_read.append('@' + bam.query_name + '\n' )
                lines_read.append(bam.query_sequence + '\n')
                lines_read.append('+' + bam.query_name + '\n')
                lines_read.append(bam.qual + '\n') """
        bamfile.close()

        bamfile = pysam.AlignmentFile(path_bam,'rb')
        #print('----------------bam length:',len_bam)
        #print('----------------total reads:',len_bam)
        
        ec_file = open(ec_dir+'/'+ec_list[k],'r')
        lines_ec = ec_file.readlines()
        tmp_ec = 0
        with open(txt_dir +'/'+ txt_list[k],'r') as lines_txt:
            time_fe = time.time()
            time_f = time_fe - time_fs  # 计算的时间差为程序的执行时间，单位为秒/s
            #print("Begin")
            time_es = time.time()
            ref_ = ''
            flag_ec = 0
            for line_txt in lines_txt:
                flag_del = list()
                line_txt_arr = line_txt.split() 
                if int(line_txt_arr[1]) == 1:
                    time_100000_b = time.time()
                if int(line_txt_arr[1]) % 100000 == 0:
                    time_100000_e = time.time()
                    #print(line_txt_arr[0],line_txt_arr[1],num_pos,time_100000_e - time_100000_b)
                    time_100000_b = time_100000_e
                    num_pos = 1
                if tmp_ec == len(lines_ec):#如果ec文件已经读取完毕，则写入剩余read
                    flag_ec = 1
                else:
                    line_ec_arr = lines_ec[tmp_ec].split()
                if ref_ != line_txt_arr[0] and ref_ != '' and flag_ec == 0:
                    if ref_ == line_ec_arr[0]:
                        #print('error type2',line_txt_arr[1],line_ec_arr[1])
                        break
                    #else:
                        #print('next ref',line_txt_arr[0])
                ref_ = line_txt_arr[0]

                if flag_ec == 0:
                    if int(line_txt_arr[1]) > int(line_ec_arr[1]) and line_txt_arr[0] == line_ec_arr[0] and flag_ec == 0:
                        #print('error type1',line_txt_arr[1],line_ec_arr[1])
                        break
                
                #构建数组
                tmp_p = 0
                bases = list(line_txt_arr[4])
                b_len = len(bases)
                l_b = 0
                for i in range(b_len):
                    if bases[l_b] == '^':
                        if len(pos) == 0:
                            pos.append([pos_tail,0])
                            pos_tail += 1
                        else :
                            #pos.append(pos[len(pos)-1]+1) #末尾元素+1
                            #pos_tail = pos[len(pos)-1] + 2
                            pos.append([pos_tail,0])
                            pos_tail += 1
                        l_b += 3
                        tmp_p += 1
                        for bam in bamfile:
                            if bam.reference_name == refs[k]:
                                bamlines.append(bam)
                                lines_read.append('@' + bam.query_name + '\n' )
                                lines_read.append(bam.query_sequence + '\n')
                                lines_read.append('+' + bam.query_name + '\n')
                                lines_read.append(bam.qual + '\n') 
                                break
                
                        #out_check.writelines(line_txt)
                    elif bases[l_b] == '$' and bases[l_b-1] != '^':
                        #tmp_p -= 1
                        #pos[tmp_p:] += 1
                        flag_del.append(tmp_p-1) #将$点位加入列表flag_del中
                        """ if tmp_p >= len(pos) or tmp_p < 0:
                            print(len(pos),tmp_p,line_txt_arr[1],count_begin,count_end) """
                        #del pos[tmp_p] 
                        l_b += 1
                        
                        #out_check.writelines(line_txt)
                    elif l_b != b_len-1 and (bases[l_b] == ',' or bases[l_b] == '.' or bases[l_b] == '*' or bases[l_b] == 'a' or bases[l_b] == 'A'or bases[l_b] == 'c' or bases[l_b] == 'C' or bases[l_b] == 'g' or bases[l_b] == 'G' or bases[l_b] == 't' or bases[l_b] == 'T' or bases[l_b] == 'N' or bases[l_b] == 'n'):
                        l_b += 1
                        tmp_p += 1
                        if bases[l_b] == '+' or bases[l_b] == '-':
                            q = 0
                            p = 1
                            l_b += 1
                            for j in range(l_b,b_len):
                                if '0'<=bases[j]<='9':
                                    l_b += 1
                                else: 
                                    break
                            for j in range(0,l_b):
                                if '0'<=bases[l_b-1]<='9' :
                                    q = q + p*int(bases[l_b-1]) + 1
                                    p = p*10
                                    l_b -= 1
                                else:   
                                    break
                            l_b = l_b + q
                        else :
                            pass
                    elif l_b == b_len-1 and (bases[l_b] == ',' or bases[l_b] == '.' or bases[l_b] == '*' or bases[l_b] == 'a' or bases[l_b] == 'A'or bases[l_b] == 'c' or bases[l_b] == 'C' or bases[l_b] == 'g' or bases[l_b] == 'G' or bases[l_b] == 't' or bases[l_b] == 'T' or bases[l_b] == 'N' or bases[l_b] == 'n'):
                        l_b += 1
                        tmp_p += 1
                    if l_b >= b_len:
                        break
                ''' while int(line_txt_arr[1]) > int(line_ec_arr[1]):
                    tmp_ec += 1
                    if tmp_ec == len(lines_ec):
                        break
                    else:
                        line_ec_arr = lines_ec[tmp_ec].split() '''

                if flag_ec == 0:
                    while line_txt_arr[1] == line_ec_arr[1] and line_txt_arr[0] == line_ec_arr[0] :
                        #print(line_ec_arr[1])
                        """ if int(line_ec_arr[2])-1 >= len(pos) or int(line_ec_arr[2])-1 < 0:
                            print(len(pos),line_ec_arr[1]) """ 
                        num = int(line_ec_arr[2])-1 #num为bam中对应read的编号,0开始
                        #bamfile = pysam.AlignmentFile(path_bam, "rb")
                        #numm = 0
                        read_name = ''
                        read_location = int(line_ec_arr[1])
                        read_ec = line_ec_arr[3]
                        read_flag = 0
                        check_location = 0
                        #print(read_ec,read_location,num,pos)
                        if read_ec == 'N' or read_ec == 'n' or line_txt_arr[2] == 'N' or line_txt_arr[2] == 'n':
                            pass
                        else:
                            read_name = bamlines[num].query_name
                            read_location = read_location - bamlines[num].pos - 1
                            read_flag = bamlines[num].flag
                            count = read_location
                            check_location = read_location
                            ci_l = 0
                            for ci in bamlines[num].cigar:
                                if ci[0] == 0: #count为比对的位置的计数 如果匹配，则-
                                    count -= ci[1]
                                elif ci[0] == 3 or ci[0] == 4:#如果有跳过 locat +
                                    read_location += ci[1]
                                elif ci[0] == 2 : #如果有删除 
                                    read_location -= ci[1]
                                    count -= ci[1]
                                elif ci[0] == 1: #如果有插入
                                    read_location += ci[1]
                                if count < 0:
                                    if line_ec_arr[4] == '*':
                                        read_location += (0-count)
                                        ci_l = ci[1] #记录碱基删除时的个数
                                    break
                            """ for bam in bamfile:
                                if numm == num:
                                    read_name = bam.query_name
                                    read_location = read_location - bam.pos - 1
                                    read_flag = bam.flag
                                    count = read_location
                                    #check_location = read_location
                                    for ci in bam.cigar:
                                        if ci[0] == 0: #count为比对的位置的计数 如果匹配，则-
                                            count -= ci[1]
                                        elif ci[0] == 3 or ci[0] == 4:#如果有跳过 locat +
                                            read_location += ci[1]
                                        elif ci[0] == 2 : #如果有删除 
                                            read_location -= ci[1]
                                            count -= ci[1]
                                        elif ci[0] == 1: #如果有插入
                                            read_location += ci[1]
                                        if count < 0:
                                            break
                                    break
                                else :
                                    numm += 1
                            bamfile.close() """
                            num_pos += 1
                            
                            l_r = 0
                            l_r = 4*num + 1 #获得该read的位置
                            str_list = list(lines_read[l_r])
                            str_list_next = list(lines_read[l_r+2])
                            ''' if read_flag == 16:
                                if read_ec == 'A' or read_ec == 'a': #比对到反链，需要改变碱基
                                    read_ec = 'T'
                                elif read_ec == 'T' or read_ec == 't':
                                    read_ec = 'A'
                                elif read_ec == 'C' or read_ec == 'c':
                                    read_ec = 'G'
                                elif read_ec == 'G' or read_ec == 'g':
                                    read_ec = 'C'
                                if len(str_list)  - 2 - read_location - pos[int(line_ec_arr[2])-1][1] < 0 or len(str_list)  - 2 - read_location - pos[int(line_ec_arr[2])-1][1] > len(str_list) - 1:
                                    print('error',read_location,len(str_list),read_name,line_ec_arr[1],num,pos) 
                                    #out_check.writelines('error' + '\t' + str(read_location) + '\t' + str(len(str_list)) + '\t' + read_name + '\t' + line_ec_arr[1] + '\t' + str(num) + '\t' + str(pos) + '\n')
                                else:
                                    if line_ec_arr[4] == '*':
                                        if(count < -5):
                                            break
                                        str_list.insert(len(str_list)  - 1 - read_location - pos[int(line_ec_arr[2])-1][1] ,read_ec)
                                        str_list_next.insert(len(str_list)  - 1 - read_location - pos[int(line_ec_arr[2])-1][1] ,'~')
                                        pos[int(line_ec_arr[2])-1][1] += 1 

                                    elif line_ec_arr[4] == '+':
                                        lenth = int(line_ec_arr[3])
                                        if (lenth > 5):
                                            break
                                        read_location += 1
                                        len_str = len(str_list)
                                        len_qual = len(str_list_next)
                                        for j in range(lenth):
                                            del str_list[len_str  - 2 - read_location - pos[int(line_ec_arr[2])-1][1] - lenth + 1]
                                            del str_list_next[len_qual  - 2 - read_location - pos[int(line_ec_arr[2])-1][1] - lenth + 1]
                                        pos[int(line_ec_arr[2])-1][1] -= lenth

                                    else:                                     
                                        str_list[len(str_list)  - 2 - read_location - pos[int(line_ec_arr[2])-1][1]] = read_ec #需要-2 因为最后是/n
                                        str_list_next[len(str_list)  - 2 - read_location - pos[int(line_ec_arr[2])-1][1]] = '~' '''
                            
                            if read_location + pos[int(line_ec_arr[2])-1][1] >= len(str_list) or read_location + pos[int(line_ec_arr[2])-1][1] < 0 :
                                print('error',read_location,check_location,len(str_list),read_name,line_ec_arr[1],num)  
                                #out_check.writelines('error' + '\t' + str(read_location) + '\t' + str(check_location) + '\t' + str(len(str_list)) + '\t' + read_name + '\t' + line_ec_arr[1] + '\t' + str(num) + '\n')
                            else:
                                if line_ec_arr[4] == '*':
                                    if(count < -5 or ci_l > 5):
                                        pass
                                    else:
                                        str_list.insert(read_location +  pos[int(line_ec_arr[2])-1][1] ,read_ec)
                                        str_list_next.insert(read_location +  pos[int(line_ec_arr[2])-1][1] ,'~')
                                        pos[int(line_ec_arr[2])-1][1] += 1  

                                elif line_ec_arr[4] == '+':
                                    lenth = int(line_ec_arr[3])
                                    if (lenth > 5):
                                        pass
                                    else:
                                        read_location += 1
                                        
                                        del str_list[read_location + pos[int(line_ec_arr[2])-1][1]:read_location + pos[int(line_ec_arr[2])-1][1]+lenth]
                                        del str_list_next[read_location + pos[int(line_ec_arr[2])-1][1]:read_location + pos[int(line_ec_arr[2])-1][1]+lenth]   
                                        pos[int(line_ec_arr[2])-1][1] -= lenth 

                                else:
                                    str_list[read_location + pos[int(line_ec_arr[2])-1][1]] = read_ec #pos[int(line_ec_arr[2])-1][1]为因add增加而调整的位置
                                    str_list_next[read_location + pos[int(line_ec_arr[2])-1][1]] = '~'
                            lines_read[l_r] = ''.join(str_list)
                            lines_read[l_r+2] = ''.join(str_list_next)
                            #lines_read[l_r][read_location] = read_ec #修改对应碱基
                        
                        tmp_ec += 1
                        if tmp_ec == len(lines_ec):
                            break
                        else:
                            line_ec_arr = lines_ec[tmp_ec].split()
            
            
                if len(flag_del) != 0:
                    for i in range(len(flag_del)):
                        del pos[flag_del[len(flag_del)-1-i]]
                        del bamlines[flag_del[len(flag_del)-1-i]]
                        out_read.writelines(lines_read[flag_del[len(flag_del)-1-i]*4:flag_del[len(flag_del)-1-i]*4+4])
                        del lines_read[flag_del[len(flag_del)-1-i]*4:flag_del[len(flag_del)-1-i]*4+4]
                    flag_del.clear()
                    
        bamfile.close()
        time_ee = time.time()
        time_e = time_ee - time_es  # 计算的时间差为程序的执行时间，单位为秒/s
        #print("ctime_e(/s): ", time_e)
        time_ws = time.time()
        #out_read.writelines(lines_read)
        time_we = time.time()
        time_w = time_we - time_ws  # 计算的时间差为程序的执行时间，单位为秒/s
        #print("ctime_w(/s): ", time_w)

    #print("Correction Finish")