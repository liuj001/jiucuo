# JiuCuo

## Overview
JiuCuo is a PacBio HiFi read correction algorithm using preassembled primary contigs based on deep image processing.

## Copy right
JiuCuo is under the Artistic License 2.0.

## Short manual

### 1. System requirements
JiuCuo is suitable for 32-bit or 64-bit machines with Linux operating systems. JiuCuo include minimap2 alignment, samtools processing and error correction, where the former two steps require relatively large memory, and the later uses no more than 10 GB.

### 2. Installation
  git clone -d master https://github.com/bbbj001/jiucuo.git
  cd jiucuo
  conda env create -f JiuCuo.yml
  pip install -r JiuCuo.txt


### 3. Inputs
- HiFi reads in FASTQ format
- Preassembled primary contigs from the reads in FASTA format

# JiuCuo

## Overview

JiuCuo is a PacBio HiFi read correction algorithm using preassembled primary contigs based on deep image processing.

## Copy right

JiuCuo is under the Artistic License 2.0.

## Short manual

### 1. System requirements

JiuCuo is suitable for 32-bit or 64-bit machines with Linux operating systems. JiuCuo include minimap2 alignment, samtools processing and error correction, where the former two steps require relatively large memory, and the later uses no more than 10 GB.

### 2. Installation
```sh
git clone -d master https://github.com/bbbj001/jiucuo.git
cd jiucuo
conda env create -f JiuCuo.yml
pip install -r JiuCuo.txt
```
### 3. Inputs
- HiFi reads in FASTQ format
- Preassembled primary contigs from the reads in FASTA format

### 4. Using JiuCuo
```sh
conda activate JiuCuo
python runJiuCuo.py -reads hifi_reads.fastq -contigs primary_contigs.fasta -output directory [-options | -options]
```
#### Mandatory:
`-reads`
  Raw HiFi reads in FASTQ format
`-contigs`
  Preassembled primary contigs from the reads in FASTA format
`-output`
  Output directory

#### Options (default value):
`-min_bases n (1)`
  Minumum number of mismatched bases required to generate an error candidate image
`-min_reads n (3)`
  Minimum number of reads required to generate the error candidate image

`-diameter_size n (100)`
  Maximum diameter size in DBSCAN
`-cluster_size n (3)`
  Minumum cluster size in DBSCAN
`-k_size n (8)`
  Size of k-mer in adapter matching
`-identity_value n (0.7)`
  Identity value in adapter matching

`-adaptor_removal (no)`
  Adapter removal from the reads
`-threads n (8)`
  Number of threads during correction
`-allocated_reads n (20000)`
  Maximum number of reads allocated to each thread
`-log (no)`
  System log to print.

### 5. Outputs
-Corrected HiFi reads

