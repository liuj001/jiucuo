# JiuCuo

## Overview

JiuCuo is a PacBio HiFi read correction algorithm using preassembled primary contigs based on deep image processing.

## Copy right

JiuCuo is under the Artistic License 2.0.

## Short manual

### 1. System requirements

JiuCuo is suitable for 32-bit or 64-bit machines with Linux operating systems. JiuCuo include minimap2 alignment, samtools processing and error correction, where the former two require relatively large memory, and the later uses no more than 10 GB.

### 2. Installation
Please ensure that Git LFS is installed. If it is not installed, please follow the steps below to install it:
```sh
sudo apt update
sudo apt install git-lfs
git lfs install
```
```sh
git clone https://github.com/liuj001/jiucuo.git
cd jiucuo
conda env create -f JiuCuo.yml
```
Notes:
- if you can't install torch and torchvision successfully, please download from this link: https://download.pytorch.org/whl/torch_stable.html
```sh
torch-1.8.1+cu111-cp38-cp38-linux_x86_64
torchvision-0.9.1+cu111-cp38-cp38-linux_x86_64
```
- The following packages may fail to install, it is recommended to install them manually.
```sh
pip install pysam
pip install lxml
pip install matplotlib
pip install pandas
pip install jsonpath
pip install scikit-learn
pip install ultralytics
pip install tqdm
conda install -c bioconda seqkit
```
### 3. Inputs
- HiFi reads in FASTQ format
- Preassembled primary contigs from the reads in FASTA format

### 4. Using JiuCuo
```sh
conda activate JiuCuo
bash runJiuCuo.sh -reads hifi_reads.fastq -contigs primary_contigs.fasta -output directory [-options | -options]
```
#### Mandatory:
`-reads`
  Raw HiFi reads in FASTQ format

`-contigs`
  Preassembled primary contigs from the reads in FASTA format

`-output`
  Output directory

#### Options (default value):
`-diameter_size [int] (600)`
  Maximum diameter size in DBSCAN

`-cluster_size [int] (3)`
  Minumum cluster size in DBSCAN

`-k_size [int] (5)`
  Size of k-mer in adapter matching

`-identity_value n (0.6)`
  Identity value in adapter matching.

`-adaptor_removal (0)`
  Adapter removal from the reads (0 is no adapter removal and 1 is adapter removal)

`-threads [int] (8)`
  Number of threads during correction

`-allocated_reads [int] (10000)`
  Maximum number of reads allocated to each thread


### 5. Outputs
`base_correction.fastq`
HiFi reads with base error correction

`base_correction_adapter_removal.fastq`
 HiFi reads with base error correction and adapter removal