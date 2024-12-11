from process import process 
from cluster_second import cluster_second
from all_detections_filt_third import base_all_detections1,ATGCFiltDetections1

import argparse
import os
import threading
from tqdm import tqdm  

def all_in_one(openfile,bam_filename,similarity):
    """
    This function runs all the processing steps sequentially on each file. 
    It executes the following steps:
    1. Processes the input file using the `process` function.
    2. Applies DBSCAN clustering using `cluster_second`.
    3. Detects adapters in sequences with `base_all_detections1`.
    4. Filters the detected adapters based on similarity using `ATGCFiltDetections1`.

    Args:
    openfile (str): The name of the file to process.
    bam_filename (str): The name of the BAM file containing sequence data.
    similarity (float): The threshold value for detecting adapter sequences based on similarity.
    """
    process(openfile) # Process the file to prepare for clustering and detection
    cluster_second(openfile,args.eps,args.min_samples) # Run DBSCAN clustering on the processed data
    base_all_detections1(openfile,bam_filename) # Detect adapter sequences in the BAM file
    ATGCFiltDetections1(openfile,similarity,args.k_size) # Filter detections based on adapter similarity
    # vote(openfile)
    

def process_names_with_threads(names, num_threads):
    """
    This function processes a list of filenames in parallel using multiple threads.
    It divides the list of files into smaller chunks and assigns each chunk to a separate thread.
    A progress bar is displayed to show the processing status.

    Args:
    names (list): A list of filenames to process.
    num_threads (int): The number of threads to use for parallel processing.
    """
    def process_names_subset(names):
        """
        This function processes a subset of filenames in a separate thread.
        
        Args:
        names (list): A subset of filenames to process.
        """
        for name in names:
            all_in_one(name,args.bam_filename,args.similarity)
            progress_bar.update(1)  # Update the progress bar after each file is processed

    threads = []
    chunk_size = len(names) // num_threads
    progress_bar = tqdm(total=len(names), desc='Processing Files')  
    # Create and start threads
    for i in range(num_threads):
        start = i * chunk_size
        end = start + chunk_size if i < num_threads - 1 else len(names)
        partial_names = names[start:end]
        t = threading.Thread(target=process_names_subset, args=(partial_names,))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    progress_bar.close() 



# Argument parser to read command-line arguments
parser = argparse.ArgumentParser(description='Process BAM and prefile files.')
parser.add_argument('--bam_filename', type=str, help='Name of the BAM file')
parser.add_argument('--similarity', type=float, help='Judge the percentage of whether it is an adapter or not. ')
parser.add_argument('--num_threads', type=int, help='Number of threads')
parser.add_argument('--eps', type=int, help='dbscan eps')
parser.add_argument('--min_samples', type=int, help='dbscan min_samples')
parser.add_argument('--k_size', type=int, help='size of k-mer')
args = parser.parse_args()
# Directory containing the files to process
directory = '../data/'

similar_filenames = []
for filename in os.listdir(directory):
    if filename.startswith('ptg') and filename.endswith('_data.csv'):
        a=str(directory)+str(filename)
        similar_filenames.append(a)
# Process the list of filenames using the specified number of threads
process_names_with_threads(similar_filenames,args.num_threads)
print("已经ok")

    
    
    
    
    