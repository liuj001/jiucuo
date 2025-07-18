import argparse
import os
import sys
import multiprocessing
from tqdm import tqdm
from process import Process
from process import process_image_name
from process import merge
from process import filter_csv
from process import images_inference
from datetime import datetime
import shutil
import glob

def get_available_threads():
    return multiprocessing.cpu_count()

def all_in_one_wrapper(args):
    
    name, bam_filename, similarity, eps, min_samples, k_size, bamview_file, adapter_output_dir = args
    try:
        
        process = Process(
            openfile=name,
            bam_filename=bam_filename,
            similarity=similarity,
            eps=eps,
            min_samples=min_samples,
            k_size=k_size,
            bamview_file=bamview_file,
            output_dir=adapter_output_dir
        )
        process.process_single_file()
        process.cluster_second()
        process.base_all_detections1()
        process.ATGCFiltDetections()
        
   
        return True
    except Exception as e:
        print(f"\nError processing {name}: {str(e)}")
        return False

def process_names_with_processes(names, args):
    
   
    task_args = [(name, args.bam_filename, args.similarity,
                 args.eps, args.min_samples, args.k_size,
                 args.bamview_file, args.adapter_output_dir) for name in names]

  
    if args.error_correction:
        st=6
    else :
        st=4
    # with tqdm(total=len(names), desc=f'STAGE {st}: Adapter removal', 
    #          position=0, bar_format="{l_bar}{bar} |") as pbar:
    for t in range(1):
        
        
        with multiprocessing.Pool(processes=args.num_threads) as pool:
          
            results = []
            for task in task_args:
                res = pool.apply_async(
                    all_in_one_wrapper,
                    args=(task,),
                # callback=lambda _: pbar.update(1)
                )
                results.append(res)
            success_count = 0
            for res in results:
                if res.get():
                    success_count += 1
            
           
            if success_count != len(names):
                print("Warning: Some files failed processing. Check error logs.")



if __name__ == "__main__":
    
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
    parser.add_argument('--error_correction', type=int, required=False, help="if run error_correction",default=1)
    
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
    
    # check min_samples is valid
    # while args.min_samples is not None and (args.min_samples < 2):
    #     print("Error: The min_samples value should more than 1, input again")
    #     args.min_samples = int(input("Please enter a (int) valid min_samples value (more than 1): "))
    
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
        # print("step 1 : start processing images")
        images_inference(args.images_dir, args.output_images_dir, args.output_detection_csv_path, args.error_correction)
        args.process_files = args.output_detection_csv_path
    else :
        args.process_files = args.output_detection_csv_path
    
    # make sure dir is existed
    os.makedirs(args.images_dir, exist_ok=True)
    os.makedirs(args.output_images_dir, exist_ok=True)
    
    
    # Directory containing the files to process
    directory = os.path.dirname(args.process_files)
    # process image name
    process_image_name(args.process_files)



    similar_filenames = []
    for filename in os.listdir(directory):
        if  filename.endswith('_data.csv'):
            a = os.path.join(directory, filename)
            similar_filenames.append(a)
    
    
    process_names_with_processes(similar_filenames, args) 
    # merge sequence adapter

    res = merge(folder_path=args.adapter_output_dir, output_name="adapter_rows_sequence.csv")
    
    # merge cluster adapter
    merge(folder_path=directory, output_name="adapter_rows_dbscan.csv", end="*second_part.csv")
    
    
    filter_csv(f'{args.adapter_output_dir}/adapter_rows_sequence.csv', f'{directory}/adapter_rows_dbscan.csv',
               f'{args.adapter_output_dir}/adapter_remove.csv')
    
    shutil.copy(f'{args.adapter_output_dir}/adapter_remove.csv', result_dir)
    
    
    for filename in os.listdir(directory):
        if filename.endswith('_data.csv'):
            abs_filepath = str(directory) + str('/') + str(filename)
            os.remove(abs_filepath)

