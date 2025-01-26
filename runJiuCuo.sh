#!/bin/bash

identity_value=0.6
threads=8
diameter_size=600
cluster_size=3
k_size=20
allocated_reads=10000
bam="/raw.bam"
bam_f="/filt.bam"
adapter_removal=0
log_file="/JIUCUO_LOG.txt"

# 参数校验函数
validate_param() {
  local value=$1
  local min=$2
  local max=$3
  local name=$4
  local is_integer=$5

  # 检查是否为整数
  if [ "$is_integer" == "1" ]; then
    if ! [[ "$value" =~ ^[0-9]+$ ]]; then
      echo "Error: option -$name needs to be an integer."
      exit 1
    fi
  else
    # 浮点数检查（仅支持简单的浮点格式）
    if ! [[ "$value" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
      echo "Error: option -$name needs to be a real number."
      exit 1
    fi
  fi

  # 使用 awk 进行小数范围比较
  if ! awk -v val="$value" -v min="$min" -v max="$max" 'BEGIN {exit !(val >= min && val <= max)}'; then
    echo "Error: option -$name needs to be between $min and $max."
    exit 1
  fi
}

# 检查 adapter_removal 是否为 0 或 1
validate_adapter_removal() {
  local value=$1
  if ! [[ "$value" =~ ^[01]$ ]]; then
    echo "Error: option -adapter_removal needs to be 0 or 1 (0 = no removal, 1 = removal)."
    exit 1
  fi
}

# 处理命令行参数
while [[ $# -gt 0 ]]; do
  case $1 in
    -contigs)     
      contigs="$2"
      shift 2 ;;
    -reads)     
      reads="$2"
      shift 2 ;;
    -output)     
      output="$2"
      shift 2 ;;
    -allocated_reads)     
      allocated_reads=$2
      shift 2 ;;
    -identity_value)     
      identity_value=$2
      shift 2 ;;
    -diameter_size)      
      diameter_size=$2
      shift 2 ;;
    -threads)      
      threads=$2
      shift 2 ;;
    -cluster_size)     
      cluster_size=$2
      shift 2 ;;
    -k_size)      
      k_size=$2
      shift 2 ;;
    -adapter_removal)      
      adapter_removal=$2
      shift 2 ;;
    *)             # 处理无效参数
      echo "echo -e "JiuCuo: PacBio HiFi read correction method using preassembled contigs based on deep image processing\nJiwen Liu, Mingfei Pan, Hongbin Wang and Ergude Bao\nGroup of Interdisciplinary Information Sciences, School of Software Engineering, Beijing Jiaotong University\n"Error: -$1 is invalid."
      exit 1 ;;
  esac
done

# 将标准输出和标准错误同时写入日志文件
exec > >(tee -a "$output$log_file") 2>&1 > /dev/null

echo -e "JiuCuo: PacBio HiFi read correction method using preassembled contigs based on deep image processing
Jiwen Liu, Mingfei Pan, Hongbin Wang and Ergude Bao
Group of Interdisciplinary Information Sciences, School of Software Engineering, Beijing Jiaotong University\n"

# 检查必要参数是否已输入
if [ -z "$contigs" ]; then
  echo "Error: mandatory -contigs is missing."
  exit 1
fi

if [ -z "$reads" ]; then
  echo "Error: mandatory -reads is missing."
  exit 1
fi

if [ -z "$output" ]; then
  echo "Error: mandatory -output is missing."
  exit 1
fi

# 参数验证
validate_param "$identity_value" 0 1 "identity_value" 0
validate_param "$diameter_size" 0 5000 "diameter_size" 1
validate_param "$cluster_size" 2 10 "cluster_size" 1
validate_param "$threads" 1 100 "threads" 1
validate_param "$k_size" 1 30 "k_size" 1
validate_param "$allocated_reads" 100 1000000 "allocated_reads" 1
validate_adapter_removal "$adapter_removal"

if [ ! -d "$output" ]; then
  mkdir "$output"
fi

echo "STAGE 1: Minimap2 alignment"
# 运行 minimap2 和 samtools
if [ ! -f "$output$bam" ]; then
  minimap2 -t32 -ax map-hifi "$contigs" "$reads" 2>> "$output/TOOLS_LOG.log" | samtools view -@ 48 -bS > "$output$bam"
fi

num=1
bamview_file="/bamview-new.csv"
bam_csv_dir="/bamview_csv"
bam_csv_dir_a="/bamview_csv/*"

python runJiuCuo.py -contigs "$contigs" -reads "$reads" -output "$output" -threads  "$threads" -allocated_reads "$allocated_reads" -adapter_removal "$adapter_removal"

if [ "$adapter_removal" -eq 0 ]; then
  echo -e "\nDone! Please find the correction file(s) in output directory."
fi

if [ "$adapter_removal" -eq "$num" ]; then
  file_count=$(ls -1 "$output$bam_csv_dir" | wc -l)
  if [ "$file_count" -eq 1 ]; then
    first_file=$(ls "$output$bam_csv_dir" | head -n 1)
    bamview_file="$bam_csv_dir/$first_file"
  else
    cat $(find "$output$bam_csv_dir" -type f -name "*.csv") > "$output$bamview_file"
  fi
fi

corr_fq="/base_correction.fastq.gz"
infile="$output$corr_fq"
corr_a_fq="/base_correction_adapter_removal.fastq.gz"
outfile="$output$corr_a_fq"
adapter="/adapter"
adapter_out="/adapter_out"
process_files="/adapter_remove.csv"
bam_m="/merge.bam"
bam_m_s="/merge_s.bam"
bam_dir="/bam"

if [ "$adapter_removal" -eq "$num" ]; then
  samtools merge "$output$bam_m" $(find "$output$bam_dir" -type f -name "*.bam") 2>> "$output/TOOLS_LOG.log"
  samtools sort "$output$bam_m" -o "$output$bam_m_s" 2>> "$output/TOOLS_LOG.log"
  samtools index "$output$bam_m_s" 2>> "$output/TOOLS_LOG.log"

  python yolo/process/run.py  --bam_filename "$output$bam_m_s" \
               --bamview_file "$output$bamview_file"\
               --images_dir "$output$adapter" \
               --similarity "$identity_value" \
               --num_threads "$threads" \
               --eps "$diameter_size" \
               --min_samples "$cluster_size"\
               --k_size "$k_size"\
               --adapter_output_dir "$output$adapter_out"

  python correction/adapter_locate-v2.py -bam "$output$bam_m_s" -outfile "$outfile" -infile "$infile" -csv "$output$adapter_out$process_files"
  echo -e "\nDone! Please find the correction file(s) in output directory."
fi
