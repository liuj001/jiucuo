# JiuCuo

## Overview

JiuCuo is a PacBio HiFi read correction algorithm using preassembled primary contigs based on deep image processing.

## Copy right

JiuCuo is under the Artistic License 2.0.

## Short manual

### 1. System requirements

JiuCuo is suitable for 32-bit or 64-bit machines with Linux Ubantu operating systems. JiuCuo usually uses no more than 40 GB memory, but needs relatively large disk space to store generated images.

### 2. Installation

Please follow the steps below to install JiuCuo.
```sh
git clone https://github.com/liuj001/jiucuo.git
pip install -U huggingface_hub
huggingface-cli download mingfeidf/JiuCuo  inceptionv4.pth  --local-dir jiucuo/inception/weight/ --local-dir-use-symlinks False --resume-download --force-download
huggingface-cli download mingfeidf/JiuCuo  best.pt  --local-dir jiucuo/yolo/process --local-dir-use-symlinks False --resume-download --force-download
cd jiucuo
conda env create -f JiuCuo.yml
conda activate JiuCuo
pip install https://download.pytorch.org/whl/cu111/torch-1.8.1%2Bcu111-cp38-cp38-linux_x86_64.whl https://download.pytorch.org/whl/cu111/torchvision-0.9.1%2Bcu111-cp38-cp38-linux_x86_64.whl -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
conda install -c bioconda seqkit
```

If the `huggingface-cli` commands fail, please run this command before them.
```sh
export HF_ENDPOINT=https://hf-mirror.com
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
`-base_correction (1) `
  Base correction in the reads (1 is base correction and 0 is no base correction)

`-adapter_removal (0)`
  Adapter removal from the reads (0 is no adapter removal and 1 is adapter removal)

`-diameter_size [int] (600)`
  Maximum diameter size in DBSCAN

`-cluster_size [int] (3)`
  Minumum cluster size in DBSCAN

`-k_size [int] (5)`
  Size of k-mer in adapter matching

`-identity_value n (0.6)`
  Identity value in adapter matching.

`-threads [int] (8)`
  Number of threads during correction

`-allocated_reads [int] (10000)`
  Maximum number of reads allocated to each thread

`-help`

### 5. Outputs
`base_correction.fastq`
HiFi reads with base error correction

`adapter_removal.fastq`
HiFi reads with adapter removal

`base_correction_adapter_removal.fastq`
 HiFi reads with base error correction and adapter removal

### 6. Fine-tuning of Inception-v4 and YOLO-v8
You could fine-tune Inception-v4 by putting your data under `inception/finetune/data`.
```sh
python inception/finetune/train.py -b 8 -lr 0.0001
```
You could fine-tune YOLO-v8 by putting your data under `yolo/finetune/data`.
```sh
pip install ultralytics==8.0.196
python yolo/finetune/train.py -b 1 -lr 0.01
```

### 7. Video demo

This [video demo](xxx) records the whole process from installation, running to obtaining the result on some tiny data.

### 8. Frequently asked questions

(1) Why is there an error during installation?

Please check and make sure the Linux system is Ubantu, the Internet is connected, and the disk usage does not exceed the limit.

(2) Why is there an error during running?

Please check and make sure the disk usage does not exceed the limit.

