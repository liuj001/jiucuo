import os 
import pandas as pd
import math

import time


def find_error(ec_dir, snp_dir, txt_dir, chr):
    time_fs = time.time()
    out_file = open(ec_dir + '/' + chr + '_ec.txt','w')
    #out_file = open('chr20_d70_ec.txt','w')
   
    sz = 0
    snp_file = snp_dir + '/' + chr + '_snp.txt'
    if os.path.exists(snp_file):
        sz = os.path.getsize(snp_file)
    if not (sz or os.path.exists(snp_file)):
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

    cname = []
    pos = []
    num = []
    base = []
    base_o = []

    with open(txt_dir + '/' + chr + '.txt','r') as f:
    #with open('chr20_d70.txt','r') as f:
        time_fe = time.time()
        time_f = time_fe - time_fs  # 计算的时间差为程序的执行时间，单位为秒/s

        print("Begin")
        time_es = time.time()
        tmp = 0
        p = 0
        for line in f:
            n_d = 0
            n = 0
            c = 0
            r = 0
            num_indel = 0
            line_arr = line.split()
            indel = 0
            j = ''
            break_ = 0

            if p == 100000:
                print(line_arr[1])
                p = 0
            else:
                p +=1

            if int(line_arr[1]) == snp[tmp] : 
                print(line_arr[1])
                if len(snp) -1 > tmp:
                    tmp += 1
                a_n = 0
                c_n = 0
                g_n = 0
                t_n = 0
                s_n = 0
                n_n = 0
                n_list = []
                m_list = []
                letter = []
                for i in line_arr[4]:
                    c +=1
                    if (i == 'a' or i == 'A') and j != '^':
                        if n_d == 0:
                            a_n +=1          
                        else:
                            n_d -=1
                    elif (i == 'c' or i == 'C') and j != '^' :
                        if n_d == 0:
                            c_n +=1          
                        else:
                            n_d -=1
                    elif (i == 'g' or i == 'G') and j != '^' :
                        if n_d == 0:
                            g_n +=1          
                        else:
                            n_d -=1
                    elif (i == 't' or i == 'T') and j != '^' :
                        if n_d == 0:
                            t_n +=1          
                        else:
                            n_d -=1
                    elif (i == 'n' or i == 'N') and j != '^' :
                        if n_d == 0:
                            n_n +=1          
                        else:
                            n_d -=1
                    elif (i == '.' or i == ',') and j != '^':
                        s_n +=1          
                    elif (i == '+' or i == '-') and j != '^':
                        indel = 1
                        num_indel += 1
                    elif indel == 1 :
                        if line_arr[4][c].isdigit():
                            n_d = (n_d+int(i))*10
                            indel = 1
                        else:
                            n_d = n_d+int(i)
                            indel = 0
                    j = i
                n_list.append(a_n)
                n_list.append(c_n)
                n_list.append(g_n)
                n_list.append(t_n)
                n_list.append(s_n)
                Max = max(n_list)
                rb_1 = 'N'
                rb_2 = 'N'
                if a_n == Max:
                    rb_1 = 'A'
                    rb_2 = 'a'
                    m_list.append(c_n)
                    m_list.append(g_n)
                    m_list.append(t_n)
                    m_list.append(s_n)
                elif c_n == Max:
                    rb_1 = 'C'
                    rb_2 = 'c'
                    m_list.append(a_n)
                    m_list.append(g_n)
                    m_list.append(t_n)
                    m_list.append(s_n)
                elif g_n == Max:
                    rb_1 = 'G'
                    rb_2 = 'g'
                    m_list.append(c_n)
                    m_list.append(a_n)
                    m_list.append(t_n)
                    m_list.append(s_n)
                elif t_n == Max:
                    rb_1 = 'T'
                    rb_2 = 't'
                    m_list.append(c_n)
                    m_list.append(g_n)
                    m_list.append(a_n)
                    m_list.append(s_n)
                elif s_n == Max:
                    rb_1 = line_arr[2]
                    m_list.append(c_n)
                    m_list.append(g_n)
                    m_list.append(t_n)
                    m_list.append(a_n)
                smax = max(m_list)
                sb_1 = 'N'
                sb_2 = 'n'
                t = int(line_arr[3])/3.0
                if smax >= t:
                    if a_n == smax:
                            sb_1 = 'A'
                            sb_2 = 'a'
                    elif c_n == smax:
                            sb_1 = 'C'
                            sb_2 = 'c'
                    elif g_n == smax:
                            sb_1 = 'G'
                            sb_2 = 'g'
                    elif t_n == smax:
                            sb_1 = 'T'
                            sb_2 = 't'
                    elif s_n == smax:
                            sb_1 = line_arr[2]
                letter.append(rb_1)
                letter.append(rb_2)
                letter.append(sb_1)
                letter.append(sb_2)
                c = 0
                for i in line_arr[4]:
                    c +=1   
                    n +=1
                    if (i == 'a' or i == 't' or i == 'c' or i == 'g' or i == 'n'or i == 'A' or i == 'T' or i == 'C' or i == 'G'or i == 'N') and j != '^':
                        if n_d == 0:
                            if i not in letter:
                                cname.append(line_arr[0])
                                pos.append(line_arr[1])
                                base_o.append(i)
                                if letter[2] == 'N':
                                    base.append(letter[0])
                                else:
                                    base.append(letter[2])
                                num.append(n)                   
                        else:
                            n_d -=1
                    elif (i == '+' or i == '-') and j != '^':
                        indel = 1
                        num_indel += 1
                        if i == '+':
                            cname.append(line_arr[0])
                            pos.append(line_arr[1])
                            base_o.append(i)
                            num.append(n-1)
                            indel = 2
                    elif indel == 1 or indel == 2:
                        #print("line_arr[4][c]=",line_arr[4][c])
                        if line_arr[4][c].isdigit():
                            n_d = (n_d+int(i))*10
                            #indel = 1
                            r +=1
                        else:
                            n_d = n_d+int(i)
                            if indel == 2:                          
                                base.append(n_d)   
                            n = n-(n_d + r+2)
                            r = 0
                            indel = 0  
                            #print("n_d=",n_d) 
                    elif (i == '.' or i == ',') and j != '^' :
                        n = n
                    elif i == '*' and j != '^':
                        cname.append(line_arr[0])
                        pos.append(line_arr[1])
                        base.append(line_arr[2])
                        base_o.append(i)
                        num.append(n)  
                    else :
                        n = n-1 
                    j = i
                    
            else:
                #print(line_arr[1])
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
                        if i == '+':
                            cname.append(line_arr[0])
                            pos.append(line_arr[1])
                            base_o.append(i)
                            num.append(n-1)
                            indel = 2
                    elif indel == 1 or indel == 2:
                        if line_arr[4][c].isdigit():
                            n_d = (n_d+int(i))*10
                            #indel = 1
                            r +=1
                        else:
                            n_d = n_d+int(i)
                            if indel == 2:                          
                                base.append(n_d)   
                            n = n-(n_d + r+2)
                            r = 0
                            indel = 0
                    elif (i == '.' or i == ',') and j != '^' :
                        n = n
                    elif i == '*' and j != '^':
                        cname.append(line_arr[0])
                        pos.append(line_arr[1])
                        base.append(line_arr[2])
                        base_o.append(i)
                        num.append(n)  
                    else :
                        n = n-1
                    j = i
    f.close()
    time_ee = time.time()
    time_e = time_ee - time_es  # 计算的时间差为程序的执行时间，单位为秒/s

    time_ws = time.time()
    for i in range(len(pos)):
        out_file.write('{:}\t{:}\t{:}\t{:}\t{:}\n'.format(cname[i],pos[i], num[i], base[i], base_o[i]))
    out_file.close()
    time_we = time.time()
    time_w = time_we - time_ws  # 计算的时间差为程序的执行时间，单位为秒/s

    print("time_e(/s): ", time_e)
    print("time_f(/s): ", time_f)
    print("time_w(/s): ", time_w)

    print("Finish")