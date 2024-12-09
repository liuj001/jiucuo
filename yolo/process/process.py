import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
import argparse
# 需要的文件有：识别所有图片中的目标的文件.csv bamview.csv
pd.options.mode.chained_assignment = None

def process(openfile):
    # 提供输入文件的路径s
    df=pd.read_csv(openfile)
    # 读入read相关信息，第一列为图片的组，第二列为read的名字，第三列为此read的其实位置
    bamview=pd.read_csv('../data/bamview-new.csv',names=['image_group','name','start'])
    # 	image_id ptg000001l-1928-68.png
    # 对于上面的这样的图片名字 插入到df里面，命名为image_ref ，取ptg000001l-1928表示这个图片的位置，是第几张图片
    df1=df['image_id'].str.rsplit('-',1).str[0] 
    # 	image_id 1	ref1078-45.png
    df.insert(7,'image_ref',df1)
    # 	bamview ref1-*	m64017_191205_225630/158075461/ccs	60583 去除重复的， 表示ptg000001l-1928 为了和上面的df对应起来
    df2=bamview['image_group'].str.rsplit('-',1).str[0].drop_duplicates()
    # 插入到bamview中 ，去除重复之后剩下的都是图片的起始位置的比对位置
    bamview.insert(3,'image_group2',df2)
    
    # 左连接 此时连接的name不是真正对应上的name
    merge_df=pd.merge(df,bamview,left_on='image_ref',right_on='image_group2',how='left')
    # 取出需要的列
    df=merge_df.iloc[:,[0,1,2,3,4,7,8,9,10]]
    
    # df.to_csv('ppp1.csv')
    # num 为 image_id的 ptg000001l-1928-68.png 中的68这个值 表示这一组read的第68张图片
    num=df['image_id'].str.split('-').str[2].str.extract(r'(\d+)').astype(int)
    num=pd.Series(num[0])
    # numY为image_id的ptg000001l-1928-68.png 中的第1928这个值 表示这一组read
    numY=df['image_id'].str.split('-').str[1].str.extract(r'(\d+)').astype(int)
    numY=pd.Series(numY[0])
    # trueX1=df['start']+(num-1)*1471+df['x1']
    trueX1=df['x1']+(num-1)*1471+df['start']
    
    trueX2=df['x2']+(num-1)*1471+df['start']
    trueY1=df['y1']+(numY-1)*901
    trueY2=df['y2']+(numY-1)*901
    
    df.insert(9,'trueX1',trueX1)
    df.insert(10,'trueY1',trueY1)
    df.insert(11,'trueX2',trueX2)
    df.insert(12,'trueY2',trueY2)
     #修正y1的高度，便于匹配是哪一条read
    # df['your_column'] = np.round(df['your_column'] / 36) * 36 变成最接近36的倍数
    fixY1=(np.round(df['y1']/36)*36).astype(int)
    num=np.round(df['y2']/36).astype(int)
    # imageNum=df['image_ref'].str[3:].astype(int)
    
    # imageNum=df['image_ref'].str.split('-',1).str[1].astype(int)
    imageNum=df['image_ref'].str.split('-').str[1].astype(int)
    imageNum=(imageNum-1)*25+num
    
    medX = (trueX1 + trueX2) / 2
    medY = (trueY1 + trueY2) / 2
    df.insert(13,'fixY1',fixY1)
    df.insert(14,'num',num)
    df.insert(15,'imageNum',imageNum)
    df.insert(16,'ImageName',imageNum)
    df.insert(17,'medX',medX)
    df.insert(18,'medY',medY)
    # df.to_csv('ppp2.csv')
    newName=np.zeros(len(df))
    # 获取read的名字 imageName read中的每一条read的名字
    specified_prefix = df['image_id'][0].split('-')[0]
        # 提取第一列中以第一个 '-' 之前的部分
        # df22=bamview.iloc[:,0]
        
    bamview['prefix'] = bamview.iloc[:,0].str.split('-').str[0]
    result_df = bamview[bamview['prefix'] == specified_prefix]
    
    for i in range(len(df)):        
        name1=result_df.iloc[df.iloc[i,15]-1,1]
        df.iloc[i,16]=name1
      
    df.to_csv(f'../data/{openfile}_p.csv')





