#补全txt文件中pos跳过的问题
import os

def add_pos(chr, txt_dir, txt_add_dir):
    out_file = open(txt_add_dir + '/' + chr + '_add.txt','w') #输出的txt文件
    len_lines = 0
    a = 0
    l = 0
    with open(txt_dir + '/' + chr + '.txt','r') as lines: #输入的txt文件
        for line in lines:
            l += 1
            line_arr = line.split()
            if int(line_arr[1]) != a + 1 and a != 0:
                print(int(line_arr[1]))
                for i in range(int(line_arr[1]) - 1 - a):
                    line2 = []
                    line2.append(line_arr[0]+ '\t'+str(a+i+1)+'   A   1   ?   !'+'\n')
                    line_ = ''.join(line2)
                    out_file.writelines(line_)
                out_file.writelines(line)
            else:
                out_file.writelines(line)
            a = int(line_arr[1])
