import os 
import pandas as pd
import math



def find_error(ec_dir, snp_dir, txt_dir, chr):
    out_file = open(ec_dir + '/' + chr + '_ec.txt','w')
   
    snp_file = snp_dir + '/' + chr + '_snp.txt'
    if os.path.exists(snp_file):
        sz = os.path.getsize(snp_file)
        if not sz:
            snp_n = []
            snp = []
            l = 0
            snp_n.append(chr)
            snp.append(l)
        else :
            ben = pd.read_csv(snp_dir + '/' + chr + '_snp.txt', sep='\t', header=None)
            snp_c = ben.iloc[: ,0]
            snp_l = ben.iloc[: ,1] 
            snp_n = list(snp_c)
            snp = list(snp_l)
    else:
        snp_n = []
        snp = []
        l = 0
        snp_n.append(chr)
        snp.append(l)

    cname = []
    pos = []
    num = []
    base = []
    base_o = []

    with open(txt_dir + '/' + chr + '.txt','r') as f:
        tmp = 0
        p = 0
        for line in f:
            n_d = 0
            n = 0
            c = 0
            r = 0
            c_d = 0
            cname_d = []
            pos_d = []
            num_d = []
            base_d = []
            base_o_d = []
            c_i = 0
            cname_i = []
            pos_i = []
            num_i = []
            base_i = []
            base_o_i = []
            num_indel = 0
            line_arr = line.split()
            indel = 0
            j = ''
            break_ = 0

            if int(line_arr[1]) == snp[tmp] : 
                if len(snp) -1 > tmp:
                    tmp += 1

            #elif int(line_arr[3]) >= 5 :        
            elif int(line_arr[3]) >= 1 :
                for i in line_arr[4]:
                    n +=1
                    c +=1
                    if (i == 'a' or i == 't' or i == 'c' or i == 'g' or i == 'n'or i == 'A' or i == 'T' or i == 'C' or i == 'G'or i == 'N') and j != '^' :
                        if n_d == 0:
                            cname.append(line_arr[0])
                            pos.append(line_arr[1])
                            base.append(line_arr[2])
                            base_o.append(i)
                            num.append(n)                   
                        else:
                            n_d -=1
                    elif (i == '+' or i == '-') and j != '^':
                        indel = 1
                        num_indel += 1
                        c_i += 1
                        if i == '+':
                            cname_i.append(line_arr[0])
                            pos_i.append(line_arr[1])
                            base_o_i.append(i)
                            num_i.append(n-1)
                            indel = 2
                    elif indel == 1 or indel == 2:
                        if line_arr[4][c].isdigit():
                            n_d = (n_d+int(i))*10
                            #indel = 1
                            r +=1
                        else:
                            n_d = n_d+int(i)
                            if indel == 2:                          
                                base_i.append(n_d)   
                            n = n-(n_d + r+2)
                            r = 0
                            indel = 0
                    elif (i == '.' or i == ',') and j != '^' :
                        n = n
                    elif i == '*' and j != '^':
                        c_d += 1                        
                        cname_d.append(line_arr[0])
                        pos_d.append(line_arr[1])
                        base_d.append(line_arr[2])
                        base_o_d.append(i)
                        num_d.append(n)  
                    else :
                        n = n-1
                    j = i
                    
                base_num = int(line_arr[3])
                #if c_d == 1 or c_d<=base_num/4:
                if c_d >= 1:
                    for i in range(len(num_d)):
                        cname.append(cname_d[i])
                        pos.append(pos_d[i])
                        base.append(base_d[i])
                        base_o.append(base_o_d[i])
                        num.append(num_d[i])   
                #if c_i == 1 or c_i<=base_num/4:
                if c_i >= 1:
                    for i in range(len(num_i)):
                        cname.append(cname_i[i])
                        pos.append(pos_i[i])
                        base.append(base_i[i])
                        base_o.append(base_o_i[i])
                        num.append(num_i[i])   
                  
    f.close()

    for i in range(len(pos)):
        out_file.write('{:}\t{:}\t{:}\t{:}\t{:}\n'.format(cname[i],pos[i], num[i], base[i], base_o[i]))
    out_file.close()
