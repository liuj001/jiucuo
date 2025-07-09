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

def adapter_pic(bam_dir, adapter_dir, chr):
    ref = chr
    
    path_pic = adapter_dir
    path_bam = bam_dir + '/' + ref + '.bam'

    # 预定义颜色映射表
    COLOR_MAP = {
        0: (255,255,255),  # 白
        1: (0,255,0),      # 绿
        2: (255,255,255),  # 白 
        3: (0,0,255),      # 蓝
        4: (255,0,0),      # 红
        5: (0,0,0),        # 黑
        7: (0,255,255)     # 青
    }
    base_template = np.zeros((901, 1471, 3), dtype=np.uint8) + 255  # 白色背景模板
    # 初始化复用数组
    img_a = np.zeros((50, 490), dtype=np.int8)
    img_A = np.zeros((901, 1471), dtype=np.int8)
    ref_n = 1
    bamlines = []
    ##读取bam
    bamfile = pysam.AlignmentFile(path_bam, "rb")
    tmp = 0
    length = 50000
    img_all = np.zeros((50, length), dtype=np.int8)
    for bam in bamfile:
        flag_save = []
        bamlines.append([bam.pos,bam.cigar])
        tmp += 1
        if tmp % 25 == 1:
            start = bam.pos
        #每25条read开始处理dui
        if tmp % 25 == 0:
            img_all.fill(0)
            #每一条处理
            for i in range(0,25):
                pos = bamlines[i][0] - start
                #每一条向后遍历
                for ci in bamlines[i][1]:
                    if pos + ci[1] >= length:
                        print('warning: read too long')
                        break
                    if ci[0] == 0:
                        img_all[i*2][pos:pos+ci[1]] = 1
                        pos += ci[1]
                    elif ci[0] == 1:
                        img_all[i*2][pos-1] = 0 - min(ci[1],45)
                        #记录生成图的位置
                        if ci[1] >= 30 and ci[1] <= 48:
                            flag_save.append(pos-1)
                    elif ci[0] == 2:
                        img_all[i*2][pos:pos+ci[1]] = 7
                        pos += ci[1]
                    # elif ci[0] == 4:
                    #     img_all[i*2][pos-1] = 0 - min(ci[1],45)
                    #     #记录生成图的位置
                    #     flag_save.append(pos-1)
                #补尾
                if bamlines[i][1][-1][0] == 4 and bamlines[i][1][-1][1] >= 30 and pos - 1 < length:
                    img_all[i*2][pos-1] = 0 - min(bamlines[i][1][-1][1],45)
                    flag_save.append(pos - 1)
                #补头
                if bamlines[i][1][0][0] == 4 and bamlines[i][1][0][1] >= 30 and bamlines[i][0] - start < length:
                    img_all[i*2][bamlines[i][0] - start] = 0 - min(bamlines[i][1][0][1],45)
                    flag_save.append(bamlines[i][0] - start)
            #判断生成图条件
            if len(flag_save) > 0:
                flag_save = sorted(list(set(flag_save)))
                j = 0
                img_n = 1
                while j < length and len(flag_save) > 0:
                    #判断是否到达生成图点位
                    if j + 444 >= flag_save[0]:
                        #先去除
                        while j + 444 >= flag_save[0]:
                            del flag_save[0]
                            if len(flag_save) == 0:
                                break
                        #再开始生成图
                        if j + 445 >= length:
                            break
                        img_a.fill(0)
                        img_a[:,:445] = img_all[:,j:j+445]
                        #如果当前没图则退出
                        if img_a.sum() == 0:
                            break
                        #将img_a中插入移到下一行
                        for k in range(0, 50, 2):  # 遍历所有奇数行（步长=2）
                            # 获取当前行和下一行的引用
                            current_row = img_a[k]
                            next_row = img_a[k+1] if k+1 < 50 else None
                            
                            if next_row is not None:
                                # 找到所有负值的列索引
                                negative_indices = np.where(current_row < 0)[0]
                                
                                # 批量计算需要设置的范围
                                lengths = -current_row[negative_indices]  # 取绝对值
                                starts = negative_indices
                                ends = np.minimum(starts + lengths, 490)  # 保证不越界
                                
                                # 使用向量化操作设置值
                                for start, end in zip(starts, ends):
                                    next_row[start:end] = 3
                                    current_row[start] = 3
    
                        #将img_a映射为img_A
                        # 第一步：处理主逻辑（img_a值非0且非2的像素）
                        i_name = str(img_n)
                        r_name = str(ref_n)
                        img_A.fill(0)
                        main_mask = (img_a != 0) & (img_a != 2)
                        i_main, j_main = np.where(main_mask)
                        
                        j1_main = 3 * j_main
                        i1_main = 18 * i_main
                        
                        # 绘制左边框
                        rows_left = i1_main[:, None] + np.arange(18)
                        rows_left = np.clip(rows_left, 0, img_A.shape[0]-1)
                        img_A[rows_left, j1_main[:, None]] = 5
                        
                        # 绘制上边框
                        cols_top = j1_main[:, None] + np.arange(3)
                        cols_top = np.clip(cols_top, 0, img_A.shape[1]-1)
                        img_A[i1_main[:, None], cols_top] = 5
                        
                        # 填充内部
                        inner_rows = i1_main[:, None] + np.arange(1,18)
                        inner_cols = j1_main[:, None] + np.arange(1,3)
                        img_A[inner_rows[:, :, None], inner_cols[:, None, :]] = img_a[i_main, j_main][:, None, None]
                        
                        # ========== 第二步：处理相邻条件 ==========
                        # 左侧相邻处理
                        valid_left_mask = (j_main >= 1)
                        left_values = img_a[i_main[valid_left_mask], j_main[valid_left_mask]-1]
                        left_cond_mask = (left_values != 0) & (left_values != 2)
                        left_neighbor_mask = np.zeros_like(j_main, dtype=bool)
                        left_neighbor_mask[valid_left_mask] = left_cond_mask
                        if left_neighbor_mask.any():
                            j_left = 3 * j_main[left_neighbor_mask]
                            i_left = 18 * i_main[left_neighbor_mask]
                            img_A[i_left[:, None] + np.arange(18), j_left[:, None]] = 5
                        
                        # 上方相邻处理
                        valid_top_mask = (i_main >= 1)
                        top_values = img_a[i_main[valid_top_mask]-1, j_main[valid_top_mask]]
                        top_cond_mask = (top_values != 0) & (top_values != 2)
                        top_neighbor_mask = np.zeros_like(i_main, dtype=bool)
                        top_neighbor_mask[valid_top_mask] = top_cond_mask
                        if top_neighbor_mask.any():
                            j_top = 3 * j_main[top_neighbor_mask]
                            i_top = 18 * i_main[top_neighbor_mask]
                            img_A[i_top[:, None], j_top[:, None] + np.arange(3)] = 5
                        
                        # ========== 第三步：常规右侧/底部边框 ==========
                        # 右侧边框
                        right_mask = (j_main < 489)
                        j_right = 3 * (j_main[right_mask] + 1)
                        i_right = 18 * i_main[right_mask]
                        img_A[i_right[:, None] + np.arange(18), j_right[:, None]] = 5
                        
                        # 底部边框
                        bottom_mask = (i_main < 49)
                        j_bottom = 3 * j_main[bottom_mask]
                        i_bottom = 18 * (i_main[bottom_mask] + 1)
                        img_A[i_bottom[:, None], j_bottom[:, None] + np.arange(3)] = 5
                        
                        # ========== 第四步：特殊边缘处理 ==========
                        # 最右侧
                        right_edge_mask = (j_main == 489)
                        if right_edge_mask.any():
                            img_A[i1_main[right_edge_mask, None] + np.arange(18), 1470] = 5
                        
                        # 最底部
                        bottom_edge_mask = (i_main == 49)
                        if bottom_edge_mask.any():
                            img_A[18*49 + 18, 3*j_main[bottom_edge_mask, None] + np.arange(3)] = 5
                        
                        # ========== 第五步：修正后的底部附加行 ==========
                        bottom_row_values = img_a[49, :]
                        valid = (bottom_row_values != 0) & (bottom_row_values != 2)
                        if valid.any():
                            j_bottom_full = 3 * np.arange(490)
                            cols_bottom_full = j_bottom_full[:, None] + np.arange(3)
                            valid_cols = cols_bottom_full[valid].flatten()
                            valid_cols = np.clip(valid_cols, 0, img_A.shape[1]-1)
                            img_A[18*49 + 18, valid_cols] = 5
                        
                        # 使用预先生成的模板加速
                        colored = base_template.copy()
                        
                        # 使用numpy布尔索引加速颜色填充
                        mask = (img_A == 1)
                        colored[mask] = (0,255,0)
                        
                        mask = (img_A == 3)
                        colored[mask] = (0,0,255)
                        
                        mask = (img_A == 4)
                        colored[mask] = (255,0,0)
                        
                        mask = (img_A == 5)
                        colored[mask] = (0,0,0)
                        
                        mask = (img_A == 7)
                        colored[mask] = (0,255,255)
        
                        # 直接保存数组
                        Image.fromarray(colored).save(path_pic + '/' + ref + '-' + r_name + '-' + i_name + '.png')
                    img_n += 1
                    j += 445
                    
            bamlines.clear()
            ref_n += 1
    
