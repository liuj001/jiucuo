import argparse
import os
import sys
import threading
from tqdm import tqdm
from process import Process
from process import process_image_name
from process import merge
from process import filter_csv
from process import images_inference
import multiprocessing
from  datetime import datetime
import shutil

def get_available_threads():
    return multiprocessing.cpu_count()



def all_in_one(openfile, bam_filename, similarity):
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
    process=Process(openfile=openfile,bam_filename=bam_filename,similarity=similarity,eps=args.eps,min_samples=args.min_samples,k_size=args.k_size,bamview_file=args.bamview_file,output_dir=args.adapter_output_dir)
    process.process_single_file().cluster_second().base_all_detections().ATGCFiltDetections()



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
            all_in_one(name, args.bam_filename, args.similarity)
            progress_bar.update(1)  # Update the progress bar after each file is processed

    threads = []
    chunk_size = len(names) // num_threads
    progress_bar = tqdm(total=len(names), desc='Processing Files',position=0)
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
parser.add_argument('--bam_filename', type=str, required=True, help='Name of the BAM file')
parser.add_argument('--bamview_file', type=str, required=False, help='Path of the bamview-new file',default="../data/bamview-new.csv")
parser.add_argument('--similarity', type=float, required=False, help='Judge the percentage of whether it is an adapter or not. ',default=0.7)
parser.add_argument('--num_threads', type=int, required=False, help='Number of threads',default=get_available_threads())
parser.add_argument('--eps', type=int, required=False,help='dbscan eps',default=100)
parser.add_argument('--min_samples', type=int, required=False, help='dbscan min_samples',default=2)
parser.add_argument('--k_size', type=int,required=False, help='size of k-mer',default=8)
parser.add_argument('--process_files', type=str,required=False, help='if not inference, setting the process file path')
parser.add_argument('--adapter_output_dir',type=str,required=False,help='adapter out put dir',default='../result')
parser.add_argument('--images_dir', type=str, required=False, help='Path of the images dir',default="../images")
parser.add_argument('--output_detection_csv_path', type=str, required=False, help='Path of  output the images inference csv file')
parser.add_argument('--output_images_dir', type=str, required=False, help='Path of the images output dir')
parser.add_argument('--is_inference', type=int, required=False, help="if inference images ,if not detection_results.csv file is needed",default=1)

args = parser.parse_args()

# get time
now = datetime.now()
data_time = "{}_{}_{}".format(now.day, now.hour, now.minute)
result_dir=args.adapter_output_dir
args.adapter_output_dir = os.path.join(args.adapter_output_dir,data_time)
os.makedirs(args.adapter_output_dir, exist_ok=True)

if args.output_images_dir is None:
    args.output_images_dir = args.adapter_output_dir
else :
    args.output_images_dir = "../output_images"


if args.output_detection_csv_path is None:
    args.output_detection_csv_path = os.path.join(args.adapter_output_dir,"detection_results.csv")
else :
    args.output_detection_csv_path = "../data/detection_results.csv"

# check min_samples is valid
while args.min_samples is not None and (args.min_samples < 2):
    print("Error: The min_samples value should more than 1, input again")
    args.min_samples = int(input("Please enter a (int) valid min_samples value (more than 1): "))

# check similarity is valid
while args.similarity is not None and (args.similarity < 0 or args.similarity > 1):
    print("Error: The similarity value should be between 0 and 1, input again")
    args.similarity = float(input("Please enter a valid similarity value (between 0 and 1): "))

# check eps is valid
while args.eps is not None and (args.eps <= 0 or args.eps >=2000 ):
    print("Error: The eps value should be between 0 and 2000, 0 is not allowed ,input again")
    args.eps = int(input("Please enter a valid eps value (between 0 and 2000 , not 0): "))

# check num_threads is valid
while args.num_threads is not None and (args.num_threads > get_available_threads()):
    print(f"cpu cores is : {get_available_threads()}, please input less than {get_available_threads()}, input again")
    args.num_threads = int(input(f"Please enter a valid number of threads (less than {get_available_threads()}): "))

# check k_size is valid
while args.k_size is None or args.k_size < 1 or args.k_size > 30:
    print("Error: The k_size value should be between 1 and 30, please input again")
    args.k_size = int(input("Please enter a valid k_size: "))

output_detection_csv_dir=os.path.dirname(args.output_detection_csv_path)
if args.is_inference > 0:
    print("step 1 : start processing images")
    images_inference(args.images_dir, args.output_images_dir, args.output_detection_csv_path)
    args.process_files = args.output_detection_csv_path

# make sure dir is existed
os.makedirs(args.images_dir, exist_ok=True)
os.makedirs(args.output_images_dir, exist_ok=True)

print("step 2 : start processing files")
# Directory containing the files to process
directory = os.path.dirname(args.process_files)
# process image name
process_image_name(args.process_files)

similar_filenames = []
for filename in os.listdir(directory):
    if  filename.endswith('_data.csv'):
        a = str(directory)+str('/') + str(filename)
        similar_filenames.append(a)
# Process the list of filenames using the specified number of threads
process_names_with_threads(similar_filenames, args.num_threads)
print("step 3 : start merge files and get results")
# merge sequence adapter
print(" start merging sequence process files ")
res = merge(folder_path=args.adapter_output_dir, output_name="adapter_rows_sequence.csv")
print(" start merging dbscan process files ")
# merge cluster adapter
merge(folder_path=directory, output_name="adapter_rows_dbscan.csv", end="*second_part.csv")


filter_csv(f'{args.adapter_output_dir}/adapter_rows_sequence.csv', f'{directory}/adapter_rows_dbscan.csv',
           f'{args.adapter_output_dir}/adapter_remove.csv')

shutil.copy(f'{args.adapter_output_dir}/adapter_remove.csv', result_dir)
print(f" results saved in {result_dir}")

for filename in os.listdir(directory):
    if filename.endswith('_data.csv'):
        abs_filepath = str(directory) + str('/') + str(filename)
        os.remove(abs_filepath)






