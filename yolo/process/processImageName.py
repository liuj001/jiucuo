import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import regex as re
import os
import csv
import pysam
import pandas as pd

# 读取CSV文件
df = pd.read_csv('../data/detection_result.csv')

# 提取第一列的图片名字
image_names = df['image_id']

# 提取第一个-之前的内容作为新的列
df['Category'] = image_names.apply(lambda x: x.split('-')[0])

# 根据Category列的唯一值创建并保存不同的CSV文件
unique_categories = df['Category'].unique()
for category in unique_categories:
    category_df = df[df['Category'] == category]

    # 去除最后一列（Category列）
    category_df = category_df.iloc[:, :-1]
    

    # 保存文件
    category_df.to_csv( f'../data/{category}_data.csv', index=False)

