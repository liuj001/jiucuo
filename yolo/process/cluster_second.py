import pandas as pd
from sklearn.cluster import DBSCAN

import numpy as np
import regex as re
from datetime import datetime
import os
import csv
import pysam
from collections import Counter

def cluster_second(openfile,eps,min_samples):
    """
    Performs clustering on image regions to identify potential adapter sequences.

    Args:
        openfile (str): The path to the processed CSV file containing image region data.

    Returns:
        None. Saves clustered noise points, which may contain adapter sequences, 
        to 'pre_adapter.csv' and 'second_part.csv'.

    Description:
        - Reads in data from '{openfile}_p.csv', focusing on median X and Y coordinates ('medX', 'medY').
        - Applies DBSCAN clustering with adjustable parameters (`eps` and `min_samples`) to detect clusters.
        - Identifies noise points (cluster label `-1`) which likely include adapters.
        - Sorts noise points by 'trueY1' for ordered alignment of adapter candidates.
        - Saves the full noise point data to '{openfile}_pre_adapter.csv'.
        - Saves selected columns of noise point data to '{openfile}_second_part.csv'.
    """
    # Read the CSV file containing image region data (contains columns 'medX', 'medY', etc.)
    df = pd.read_csv(f'{openfile}_p.csv')
    
    # Extract 'medX' and 'medY' features from the DataFrame for clustering
    features = df[['medX','medY']]
    
    # Apply DBSCAN clustering with the provided parameters eps and min_samples
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    labels = dbscan.fit_predict(features)
    
    # Identify the noise points (label == -1) which are likely to contain adapter sequences
    noise_points_original = df.loc[labels == -1].copy()
    
    # Identify the points that are part of a cluster (not noise, label != -1)
    cluster=df.loc[~(labels==-1)].copy()
    
    # Sort the noise points and cluster points by the 'trueY1' value for proper ordering
    cluster.sort_values(by='trueY1',ascending=True)
    noise_points_original.sort_values(by='trueY1', ascending=True)
    
    #  Save all noise points (potential adapter sequences) to a CSV file
    noise_points_original.to_csv(f'{openfile}_pre_adapter.csv')
    
    #Add additional columns to the cluster DataFrame for future processing
    cluster['confidence']=0
    cluster['hard_score']=0
    cluster_process=cluster[['image_id','x1','y1','x2','y2','image_ref','image_group','name','start',
                             'trueX1','trueY1','trueX2','trueY2','fixY1','num','imageNum','ImageName','medX','medY']]
    # Save the cluster data to a CSV file for further investigation
    cluster_process.to_csv(f'{openfile}_cluster.csv')
    noise_points_original_columns_to_save=noise_points_original.columns[1:]
    data_to_save = noise_points_original[noise_points_original_columns_to_save]
    openfile=openfile[8:-4]
    #Save the selected noise point data to a CSV file for further analysis
    data_to_save.to_csv(f'../data/{openfile}_second_part.csv',index=False)
    
