import glob
import os
import numpy as np
import pysam
import pandas as pd
import csv
from sklearn.cluster import DBSCAN
from ultralytics import YOLO
import cv2
from tqdm import tqdm


def images_inference(image_dir, output_image_dir, output_csv_path):
    # Load the pre-trained YOLO model
    model = YOLO("yolo/process/best.pt")
    image_paths = [os.path.join(image_dir, fname) for fname in os.listdir(image_dir) if
                   fname.endswith(('.jpg', '.png', '.jpeg'))]

    # Output directory to save results
    os.makedirs(output_image_dir, exist_ok=True)

    # Store detection results
    all_results = []

    # Use tqdm to show a progress bar for the loop
    for image_path in tqdm(image_paths, desc="Processing images", unit="image"):  # Add tqdm progress bar
        result = model(image_path, verbose=False)  # Run inference on a single image
        image_name = os.path.basename(image_path)
        boxes = result[0].boxes
        if boxes is not None:
            for box in boxes:
                # Convert Tensor data to standard Python types
                x1, y1, x2, y2 = map(int, box.xyxy[0].tolist())
                confidence = float(box.conf[0])
                class_id = int(box.cls[0])
                class_name = model.names[class_id]
                all_results.append([image_name, x1, y1, x2, y2, confidence, class_id, class_name])

        # Save images
        # output_image_path = os.path.join(output_image_dir, image_name)
        # annotated_image = result[0].plot()
        # cv2.imwrite(output_image_path, annotated_image)

    # Save detection results to a CSV file
    df = pd.DataFrame(all_results,
                      columns=["image_id", "x1", "y1", "x2", "y2", "confidence", "hard_score", "Class Name"])
    new_df = df.iloc[:, :7]
    new_df.to_csv(output_csv_path, index=False, encoding="utf-8-sig")

    print(f"Detection results saved to {output_csv_path}")
    # print(f"Annotated images saved in {output_image_dir}")


def process_image_name(process_image_name_file="../data/detection_results.csv"):
    # remove ptg* files
    process_image_name_dir = os.path.dirname(process_image_name_file)
    for filename in os.listdir(process_image_name_dir):
        if filename.endswith('_data.csv'):
            a = str(process_image_name_dir) + str('/') + str(filename)
            os.remove(a)

    # for file in glob.glob(os.path.join(process_image_name_dir, "chr*")):
    #     os.remove(file)
    df = pd.read_csv(process_image_name_file)
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
        # category_df.to_csv(f'../data/{category}_data.csv', index=False)
        category_df.to_csv(f'{process_image_name_dir}/{category}_data.csv', index=False)


def filter_csv(input_file_1, input_file_2, output_file):
    # read the first file
    a = set()
    with open(input_file_1, mode='r', newline='', encoding='utf-8') as f1:
        reader = csv.reader(f1)
        for row in reader:
            a.add(f"{row[0]},{row[1]}")

    outname = output_file.rsplit(".", 1)
    adapter_out_file = outname[0] + str('_read_name.csv')

    # read the second file and filter
    with open(input_file_2, mode='r', newline='', encoding='utf-8') as f2, \
            open(output_file, mode='w', newline='', encoding='utf-8') as f_out, \
            open(adapter_out_file, mode='w', newline='', encoding='utf-8') as f_adapter_name:
        reader = csv.reader(f2)
        writer = csv.writer(f_out)
        writer1 = csv.writer(f_adapter_name)

        for row in reader:
            # combine the first file first key and the second key
            if f"{row[0]},{row[1]}" in a:
                writer.writerow(row[:5])
                writer1.writerow([row[16]])

    print(f"Result data has been saved to {output_file}")


def merge(folder_path="../result", output_name="adapter_rows_sequence.csv", end="*all_detections_to_pre_adapter.csv"):
    output_file = os.path.join(folder_path, output_name)

    header_written = False
    if not glob.glob(os.path.join(folder_path, end)):

        header = ['image_id', 'x1', 'y1', 'x2', 'y2']
        with open(output_file, 'w', newline='') as file:
            writer = csv.writer(file)
            # 写入表头
            writer.writerow(header)

    for file in glob.glob(os.path.join(folder_path, end)):

        if os.path.exists(file):
            with open(file, "r") as f_in:
                lines = f_in.readlines()

                if not header_written:
                    with open(output_file, "w") as f_out:
                        f_out.write(lines[0])
                    header_written = True

                with open(output_file, "a") as f_out:
                    f_out.writelines(lines[1:])
    print(f"Saved {output_name} in : {output_file}")
    # 删除所有 *all_detections_to_pre_adapter.csv 文件
    for file in glob.glob(os.path.join(folder_path, end)):
        os.remove(file)


class Process:
    def __init__(self, openfile, bam_filename, similarity, eps, min_samples, k_size,
                 bamview_file="../data/bamview-new.csv", output_dir='../result'):
        # Disable pandas warning
        pd.options.mode.chained_assignment = None
        # init args
        self.openfile = openfile
        self.data_name = os.path.basename(openfile)
        self.dir_name = os.path.dirname(openfile)
        self.bamview_new_file = bamview_file
        self.bam_filename = bam_filename
        self.similarity = similarity
        self.eps = eps
        self.min_samples = min_samples
        self.k_size = k_size
        self.output_dir = output_dir

    def process_single_file(self):
        # Provide the input file path
        df = pd.read_csv(self.openfile)
        # Read related information for reads, where the first column is the image group,
        # the second column is the read name, and the third column is the start position of the read
        bamview = pd.read_csv(self.bamview_new_file, names=['image_group', 'name', 'start'])

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
        df.sort_values(by='trueY1', ascending=True, inplace=True)
        self.df = df
        return self

    def cluster_second(self):
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
        # df = pd.read_csv(f'{self.openfile}_p.csv')

        # Extract 'medX' and 'medY' features from the DataFrame for clustering
        df = self.df
        features = df[['medX', 'medY']]

        # Apply DBSCAN clustering with the provided parameters eps and min_samples
        dbscan = DBSCAN(eps=self.eps, min_samples=self.min_samples)
        labels = dbscan.fit_predict(features)

        # Identify the noise points (label == -1) which are likely to contain adapter sequences
        noise_points_original = df.loc[labels == -1].copy()

        noise_points_original_columns_to_save = noise_points_original.columns[:]
        data_to_save = noise_points_original[noise_points_original_columns_to_save]

        # Save the selected noise point data to a CSV file for further analysis
        data_to_save.to_csv(f'{self.dir_name}/{self.data_name}_second_part.csv', index=False)
        return self

    def find_substring_indices(self, long_string, string_array, window_size):
        """
        This function finds all occurrences of substrings in a given long string.
        The substrings should be of the specified window_size and present in string_array.

        Args:
        long_string (str): The string in which we want to find the substrings.
        string_array (list): List of substrings to search for.
        window_size (int): The length of each substring to find.

        Returns:
        dict: A dictionary where the key is the substring, and the value is a list of indices where the substring is found.
        """
        substring_indices = {}
        substring_length = window_size
        for i in range(len(long_string) - substring_length + 1):
            substring = long_string[i:i + substring_length]
            if substring in string_array:
                if substring not in substring_indices:
                    substring_indices[substring] = []
                substring_indices[substring].append(i)
        return substring_indices

    def find_substring_indices(self, long_string, string_array, window_size):
        """
        This function finds all occurrences of substrings in a given long string.
        The substrings should be of the specified window_size and present in string_array.

        Args:
        long_string (str): The string in which we want to find the substrings.
        string_array (list): List of substrings to search for.
        window_size (int): The length of each substring to find.

        Returns:
        dict: A dictionary where the key is the substring, and the value is a list of indices where the substring is found.
        """
        substring_indices = []
        substring_length = window_size
        for i in range(len(long_string) - substring_length + 1):
            indices = []
            substring = long_string[i:i + substring_length]
            for j in range(len(string_array)):
                if string_array[j] == substring:
                    indices.append(j)
            if len(indices) > 0:
                substring_indices.append(indices)
        return substring_indices

    def process_dict(self, input_dict):
        """
        This function processes a dictionary of lists of indices.
        It merges consecutive indices and returns a flattened list of non-consecutive indices.

        Args:
        input_dict (dict): Dictionary with lists of indices.

        Returns:
        list: A list of processed indices.
        """
        # Convert the values of the dictionary into an array list
        array_list = input_dict

        # Initialize the result array
        result = []

        for i in range(len(array_list)):

            current_array = array_list[i]
            if i == 0:

                # For the first array, directly append the first element
                result.append(current_array[0])
            else:
                previous = result[-1]
                found = False
                # Merge consecutive numbers
                for num in current_array:
                    if num == previous + 1:
                        result.append(num)
                        found = True
                        break
                # If no consecutive numbers found, select the smallest number greater than the previous one
                if not found:
                    for num in current_array:
                        if num > previous:
                            result.append(num)
                            break

        return result

    def calculate_result(self, arr):
        """
        This function calculates a result value based on the distances between consecutive elements in the list `arr`.

        Args:
        window_size (int): The size of the window to calculate over.
        arr (list): A list of indices or distances to process.

        Returns:
        int: The calculated result based on the conditions.
        """
        result = 0
        n = len(arr)
        for i in range(n):
            if i == 0:
                result += self.k_size
            elif arr[i] - arr[i - 1] >= self.k_size:
                result += self.k_size
            elif arr[i] - arr[i - 1] < self.k_size:
                x = arr[i] - arr[i - 1]
                result += x
        return result

    def union_func(self, long_string, string_array, len):
        """
        This function calculates a union result based on substring matches, using previously defined helper functions.

        Args:
        long_string (str): The long string to check substrings against.
        string_array (list): The list of substrings to check for.
        window_size (int): The size of the window to match substrings.
        len (int): The length used for normalizing the result.

        Returns:
        float: A rounded value of the calculated result normalized by `len`.
        """
        indices = self.find_substring_indices(long_string, string_array, self.k_size)
        arr = self.process_dict(indices)
        calculated_result = self.calculate_result(arr)
        return round(calculated_result / len, 2)

    def contains_adapter_check(self, x, thresholds, new_array1, new_array2, new_array3, new_array4):
        """
        This function checks if the input sequence contains any adapter sequence based on pre-defined threshold values.

        Args:
        x (str): The sequence to check.
        window_size (int): The size of the window to check.
        thresholds (dict): A dictionary of thresholds for the various adapter arrays.
        new_array1, new_array2, new_array3, new_array4 (list): The pre-processed adapter sequences.

        Returns:
        bool: Returns True if an adapter sequence is found; otherwise, False.
        """

        for new_array, threshold in thresholds.items():
            count = 0
            if len(x) == 35:
                if new_array == 'new_array2':
                    count = self.union_func(x, new_array2, len=35)
                if new_array == 'new_array4':
                    count = self.union_func(x, new_array4, len=35)
            if len(x) == 45:
                if new_array == 'new_array1':
                    count = self.union_func(x, new_array1, len=45)
                if new_array == 'new_array3':
                    count = self.union_func(x, new_array3, len=45)

            if count >= threshold:
                return True

        return False

    def base_all_detections(self):
        """
        This function processes the BAM file and matches the sequences from the CSV file to identify reads containing adapter sequences.

        Args:
        openfile (str): The input file to process.
        bam_filename (str): The BAM file containing the sequence reads.

        Returns:
        pd.DataFrame: A DataFrame containing rows with detected adapter sequences.
        """
        out_file = open(f'{self.data_name}_pre_adapter_out_all_detections.csv', 'w')

        self.df.to_csv(f'{self.data_name}_tmp.csv', index=False)
        with open(f'{self.data_name}_tmp.csv', 'r') as in_file:
            reader = csv.reader(in_file)
            rows = list(reader)
        writer = csv.writer(out_file)
        rows[0].append('read_seq')
        writer.writerow(rows[0])
        tmp = 1
        bamfile = pysam.AlignmentFile(self.bam_filename, "rb")

        for bam in bamfile:
            """ if bam.pos == 60586:
                print(bam.cigar)
                break """
            if len(rows[tmp]) < 16:
                print(rows[tmp])

            while rows[tmp][16] == bam.query_name:
                row = rows[tmp]
                start = int(row[8])
                r2 = row[0].split('-')
                r22 = r2[2].split('.')

                num = int(r22[0]) - 1
                truex1 = start + num * 445 + int((int(row[1]) - 1) / 3)
                truex2 = start + num * 445 + int((int(row[3]) - 1) / 3)
                count = 0
                locat = 0
                for ci in bam.cigar:
                    if ci[0] == 4 or ci[0] == 1:
                        count += ci[1]
                    elif ci[0] == 2:
                        count -= ci[1]
                        locat += ci[1]
                    elif ci[0] == 0:
                        locat += ci[1]
                    if locat >= truex1 - bam.pos:
                        break
                if truex1 - bam.pos <= 1:  # Head of the sequence
                    count -= 45
                elif truex2 - bam.pos >= locat + 30:  # Tail of the sequence
                    count += 1
                truex1 += count
                truex2 += count
                if truex1 - bam.pos < 0:
                    add = bam.pos - truex1
                    truex1 += add
                    truex2 += add
                if truex1 - bam.pos - 8 >= 0:
                    left = truex1 - bam.pos
                else:
                    left = truex1 - bam.pos
                sequence = bam.query_sequence[left:truex2 - bam.pos + 8]

                row.append(sequence)
                writer = csv.writer(out_file)
                writer.writerow(row)
                tmp += 1
                if tmp == len(rows):
                    break
            if tmp == len(rows):
                break
        for file in glob.glob(f'{self.data_name}_tmp.csv'):
            os.remove(file)

        return self

    def ATGCFiltDetections(self):
        """
        This function performs detection on adapter sequences based on similarity and k-mer size,
        filtering sequences based on the presence of adapter sequences.

        Args:
        openfile (str): The input file to process.
        similarity (float): The threshold similarity value for detecting adapter sequences.
        k_size (int): The window size for matching adapter sequences.

        Returns:
        pd.DataFrame: A DataFrame containing rows with detected adapter sequences.
        """
        # Define the adapter sequences
        self.sequence_df = pd.read_csv(f'{self.data_name}_pre_adapter_out_all_detections.csv')
        df = self.sequence_df
        adapter_sequences = ['ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT', 'AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA',
                             'ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT', 'TCCTCCTCCTCCGTTAATTTTTTTTTTTTTTTTTT']
        # Prepare the arrays for each adapter sequence
        new_array1 = []
        new_array2 = []
        new_array3 = []
        new_array4 = []
        window = self.k_size
        sequence0 = adapter_sequences[0]
        for i in range(len(sequence0) - window + 1):
            new_array1.append(sequence0[i:i + window])

        sequence1 = adapter_sequences[1]
        for i in range(len(sequence1) - window + 1):
            new_array2.append(sequence1[i:i + window])

        sequence2 = adapter_sequences[2]
        for i in range(len(sequence2) - window + 1):
            new_array3.append(sequence2[i:i + window])

        sequence3 = adapter_sequences[3]
        for i in range(len(sequence3) - window + 1):
            new_array4.append(sequence3[i:i + window])
        # Apply the adapter check for each sequence
        df['read_seq'] = df['read_seq'].astype(str)

        thresholds = {'new_array1': self.similarity, 'new_array2': self.similarity, 'new_array3': self.similarity,
                      'new_array4': self.similarity}
        window_size = self.k_size
        size = 45
        df['contains_adapter'] = df['read_seq'].apply(
            lambda seq: any(
                self.contains_adapter_check(sub_seq, thresholds, new_array1, new_array2, new_array3, new_array4)
                for sub_seq in [seq[i:i + size] for i in range(len(seq) - size + 1)]
            )
        )

        size1 = 35
        df['contains_adapter1'] = df['read_seq'].apply(
            lambda seq: any(
                self.contains_adapter_check(sub_seq, thresholds, new_array1, new_array2, new_array3, new_array4)
                for sub_seq in [seq[i:i + size1] for i in range(len(seq) - size1 + 1)]
            )
        )
        # Filter and save rows containing adapter sequences

        adapter_rows = df[df['contains_adapter'] == True]
        adapter_rows1 = df[df['contains_adapter1'] == True]

        file_name = f'all_detections_to_pre_adapter.csv'
        adapter_rows_columns_to_save = adapter_rows.columns[0:]
        data_to_save = adapter_rows[adapter_rows_columns_to_save]

        adapter_rows1_columns_to_save = adapter_rows1.columns[0:]
        data_to_save1 = adapter_rows1[adapter_rows1_columns_to_save]

        combined_data = pd.concat([data_to_save1, data_to_save])
        # print(combined_data)
        combined_data = combined_data.drop_duplicates(subset=['image_id', 'x1'])
        # Save the final result
        # combined_data.to_csv(f'{self.openfile}_{file_name}', index=False)

        if not combined_data.empty:
            combined_data.to_csv(f'{self.output_dir}/{self.data_name}_{file_name}', index=False)

        os.remove(f'{self.data_name}_pre_adapter_out_all_detections.csv')
        return self
