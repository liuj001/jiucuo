import os #生成候选文件

def candidate_make(chr,txt_dir, bcf_txt_dir, min_bases, min_reads):
    out_file2 = open(bcf_txt_dir+'/'+chr+'.bcf.txt','w') #候选txt文件
    with open(txt_dir+'/'+chr+'.txt','r') as f: #原txt文件
        for line in f:
            num = 0
            num_str = ''
            line_arr = line.split()
            base_num = int(line_arr[3])
            indel = 0
            j = ''
            if line_arr[2] != 'N'  :
                for i in line_arr[4]:
                    if i == 'a' or i == 't' or i == 'c' or i == 'g' or i == 'A' or i == 'T' or i == 'C' or i == 'G' :
                        num += 1
                        indel = 0
                        if num_str != '':
                            num -= int(num_str)
                            num_str = ''
                    elif (i == '+' or i == '-') and j != '^':
                        indel = 1
                        
                    elif indel == 1 :
                        num_str += i
                    elif i == '*' and j != '^':
                        base_num -= 1
                    j = i
                #if (base_num <= 3 and num >= 2) or (base_num > 3 and num >= base_num/3 ) : #2024/04/18改动
                if (base_num >= min_reads and num >= min_bases):
                #if (base_num >= 10 and num >= base_num/3 ):
                #if (base_num >= 15 and num >= base_num/3 ):
                    line1 = ''.join(line_arr[0])
                    line2 = ''.join(line_arr[1])
                    out_file2.writelines(line1  + '\t' +  line2 + '\n' )
                    