import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
import argparse

# Disable pandas warning
pd.options.mode.chained_assignment = None

def process(openfile):
    # Provide the input file path
    df = pd.read_csv(openfile)
    # Read related information for reads, where the first column is the image group, 
    # the second column is the read name, and the third column is the start position of the read
    bamview = pd.read_csv('../data/bamview-new.csv', names=['image_group', 'name', 'start'])
    
    # Process the image names, inserting into df and creating 'image_ref', which extracts the image number 
    # (e.g., ptg000001l-1928-68.png -> ptg000001l-1928)
    df1 = df['image_id'].str.rsplit('-', n=1).str[0] 
    df.insert(7, 'image_ref', df1)
    
    # Process bamview image names to remove duplicates and get the starting alignment position of images
    df2 = bamview['image_group'].str.rsplit('-', n=1).str[0].drop_duplicates()
    bamview.insert(3, 'image_group2', df2)
    
    # Perform a left join to merge image-related information with bamview data
    merge_df = pd.merge(df, bamview, left_on='image_ref', right_on='image_group2', how='left')
    
    # Extract the necessary columns
    df = merge_df.iloc[:, [0, 1, 2, 3, 4, 7, 8, 9, 10]]
    
    # 'num' represents the image number extracted from the image ID, indicating the nth image in the read group
    num = df['image_id'].str.split('-').str[2].str.extract(r'(\d+)').astype(int)
    num = pd.Series(num[0])
    
    # 'numY' represents the read group extracted from the image ID
    numY = df['image_id'].str.split('-').str[1].str.extract(r'(\d+)').astype(int)
    numY = pd.Series(numY[0])
    
    # Calculate true coordinates: trueX1, trueX2, trueY1, trueY2
    trueX1 = df['x1'] + (num - 1) * 1471 + df['start']
    trueX2 = df['x2'] + (num - 1) * 1471 + df['start']
    trueY1 = df['y1'] + (numY - 1) * 901
    trueY2 = df['y2'] + (numY - 1) * 901
    
    # Insert true coordinates into df
    df.insert(9, 'trueX1', trueX1)
    df.insert(10, 'trueY1', trueY1)
    df.insert(11, 'trueX2', trueX2)
    df.insert(12, 'trueY2', trueY2)
    
    # Fix y1 height to the nearest multiple of 36
    fixY1 = (np.round(df['y1'] / 36) * 36).astype(int)
    
    # Calculate and assign the 'num' and 'imageNum' columns
    num = np.round(df['y2'] / 36).astype(int)
    imageNum = df['image_ref'].str.split('-').str[1].astype(int)
    imageNum = (imageNum - 1) * 25 + num
    
    # Calculate the median values for X and Y coordinates
    medX = (trueX1 + trueX2) / 2
    medY = (trueY1 + trueY2) / 2
    
    # Insert calculated columns into df
    df.insert(13, 'fixY1', fixY1)
    df.insert(14, 'num', num)
    df.insert(15, 'imageNum', imageNum)
    df.insert(16, 'ImageName', imageNum)
    df.insert(17, 'medX', medX)
    df.insert(18, 'medY', medY)
    
    # Initialize a new column for read names
    newName = np.zeros(len(df))
    
    # Get the read name for each read based on its image group
    specified_prefix = df['image_id'][0].split('-')[0]
    
    # Add prefix to bamview data for easier matching
    bamview['prefix'] = bamview.iloc[:, 0].str.split('-').str[0]
    result_df = bamview[bamview['prefix'] == specified_prefix]
    
    # For each row in df, assign the corresponding read name
    for i in range(len(df)):        
        name1 = result_df.iloc[df.iloc[i, 15] - 1, 1]
        df.iloc[i, 16] = name1
    
    # Save the processed dataframe to a new CSV file
    df.to_csv(f'../data/{openfile}_p.csv')
