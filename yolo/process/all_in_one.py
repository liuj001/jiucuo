from process import process 
from cluster_second import cluster_second
from all_detections_filt_third import base_all_detections1,ATGCFiltDetections1

import argparse
import os
import threading
from tqdm import tqdm  

def all_in_one(openfile,bam_filename,similarity):
    process(openfile)
    cluster_second(openfile,args.eps,args.min_samples)
    base_all_detections1(openfile,bam_filename)
    ATGCFiltDetections1(openfile,similarity,args.k_size)
    # vote(openfile)
    

def process_names_with_threads(names, num_threads):
    def process_names_subset(names):
        for name in names:
            all_in_one(name,args.bam_filename,args.similarity)
            progress_bar.update(1)  

    threads = []
    chunk_size = len(names) // num_threads
    progress_bar = tqdm(total=len(names), desc='Processing Files')  
    
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




parser = argparse.ArgumentParser(description='Process BAM and prefile files.')
parser.add_argument('--bam_filename', type=str, help='Name of the BAM file')
parser.add_argument('--similarity', type=float, help='Judge the percentage of whether it is an adapter or not. ')
parser.add_argument('--num_threads', type=int, help='Number of threads')
parser.add_argument('--eps', type=int, help='dbscan eps')
parser.add_argument('--min_samples', type=int, help='dbscan min_samples')
parser.add_argument('--k_size', type=int, help='size of k-mer')
args = parser.parse_args()

directory = '../data/'

similar_filenames = []
for filename in os.listdir(directory):
    if filename.startswith('ptg') and filename.endswith('_data.csv'):
        a=str(directory)+str(filename)
        similar_filenames.append(a)

process_names_with_threads(similar_filenames,args.num_threads)

    
    
    
    
    