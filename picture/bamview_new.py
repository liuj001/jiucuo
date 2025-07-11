import csv
import os

def bamview(csv_f,txt_f,bamview_txt_dir):
    #open_file = open('bamview.txt','r')
    
    out_file = open(csv_f, 'w')
    txt=bamview_txt_dir+'/'+txt_f
    #lines = open_file.readlines()
    a = 0
    b = 1
    str_ = ''
    with open(txt,'r') as lines:
        for line in lines:
            if line[0] == '@':
                pass
            else :
                line_arr = line.split()
                if str_ != line_arr[2] and str_ != '':
                    b = 1
                    a = 0
                str_ = line_arr[2]
                writer = csv.writer(out_file)
                data = [txt_f.replace('_bamview.txt','') + '-' + str(b) + '-*',line_arr[0],line_arr[3]]
                writer.writerow(data)
                a += 1
                if a == 25:
                    a = 0
                    b += 1

    out_file.close()