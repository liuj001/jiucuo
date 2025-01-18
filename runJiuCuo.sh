#!/bin/bash
identity_value=0.6
threads=8
diameter_size=600
cluster_size=3
k_size=20
#min_bases=1
#min_reads=3
allocated_reads=10000
bam="/raw.bam"
bam_f="/filt.bam"
adaptor_removal=0

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
      echo "Error: $name must be an integer."
      exit 1
    fi
  else
    # 浮点数检查（仅支持简单的浮点格式）
    if ! [[ "$value" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
      echo "Error: $name must be a valid number."
      exit 1
    fi
  fi

# 使用 awk 进行小数范围比较
  if ! awk -v val="$value" -v min="$min" -v max="$max" 'BEGIN {exit !(val >= min && val <= max)}'; then
    echo "Error: $name must be between $min and $max."
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
    -adaptor_removal)      
      adaptor_removal=1
      shift 2 ;;
    *)             # 处理无效参数
      echo "Invalid option: $1"
      exit 1 ;;
  esac
done

validate_param "$identity_value" 0 1 "identity_value" 0
validate_param "$diameter_size" 0 5000 "diameter_size" 0
validate_param "$cluster_size" 2 10 "cluster_size" 1
validate_param "$threads" 1 100 "threads" 1
validate_param "$k_size" 1 30 "k_size" 1
validate_param "$allocated_reads" 100 1000000 "allocated_reads" 1

if [ ! -d "$output" ]; then
  mkdir "$output"
fi

# 运行 minimap2 和 samtools
if [ ! -f "$output$bam" ]; then
  #minimap2 -t32 -ax map-hifi "$contigs" "$reads" | samtools view -@ 48 -bS > "$output$bam"
  minimap2 -t32 -ax map-hifi "$contigs" "$reads" 2>> "$output/runJiuCuo.log" | samtools view -@ 48 -bS 2>> "$output/runJiuCuo.log" > "$output$bam"

fi

num=1
bam_txt="/bamview.txt"
bamview_file="/bamview-new.csv"

if [ "$adaptor_removal" -eq "$num" ]; then
  samtools view $output$bam_f 2>> "$output/runJiuCuo.log" > "$output$bam_txt"
  python picture/bamview-new.py -cvs "$output$bamview_file" -txt "$output$bam_txt"
fi

python runJiuCuo.py -contigs "$contigs" -reads "$reads" -output "$output" -threads  "$threads" -min_bases "$min_bases" -min_reads "$min_reads" -allocated_reads "$allocated_reads" -adaptor_removal "$adaptor_removal"

corr_fq="/correction.fastq.gz"
infile="$output$corr_fq"
corr_a_fq="/correction_ar.fastq.gz"
outfile="$output$corr_a_fq"
adapter="/adapter"
adapter_out="/adapter_out"
process_files="/adapter_remove.csv"

if [ "$adaptor_removal" -eq "$num" ]; then
  python yolo/process/run.py  --bam_filename "$output$bam_f" \
               --bamview_file "$output$bamview_file"\
               --images_dir "$output$adapter" \
               --similarity "$identity_value" \
               --num_threads "$threads" \
               --eps "$diameter_size" \
               --min_samples "$cluster_size"\
               --k_size "$k_size"\
               --adapter_output_dir "$output$adapter_out"\
  python correction/adapter_locate-v2.py -bam "$output$bam_f" -outfile "$outfile" -infile "$infile" -csv "$output$adapter_out$process_files"
fi
