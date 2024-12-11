import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import regex as re
import os
import csv
import pysam
import pandas as pd

# Read the CSV file
df = pd.read_csv('../data/detection_result.csv')

# Extract the image names from the first column
image_names = df['image_id']

# Extract the part before the first '-' and create a new column 'Category'
df['Category'] = image_names.apply(lambda x: x.split('-')[0])

# Create and save separate CSV files based on the unique values in the 'Category' column
unique_categories = df['Category'].unique()
for category in unique_categories:
    category_df = df[df['Category'] == category]

    # Remove the last column ('Category' column)
    category_df = category_df.iloc[:, :-1]
    
    # Save the data to a new CSV file
    category_df.to_csv(f'../data/{category}_data.csv', index=False)
