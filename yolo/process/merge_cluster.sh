#!/bin/bash
# 该文件用于合并所有在../result目录下contig的输出结果，输出为adapter_rows.csv
cd ../data
# 指定包含CSV文件的文件夹路径
folder_path="."

# 创建一个空的CSV文件
output_file="$folder_path/adapter_rows_second.csv"
touch "$output_file"

# 获取当前目录下以 "adapter_rows" 开头的第一个CSV文件
first_file=$(find "$folder_path" -maxdepth 1 -type f -name '*second_part.csv' | sort | head -n 1)
echo "$first_file"
cat $first_file
# 如果找到第一个文件，则将其表头添加到输出文件中
if [ -n "$first_file" ]; then
    echo "--------------------"
    head -n 1 "$first_file" > "$output_file"
    # cat $output_file
fi

# 循环读取以 "adapter_rows" 开头的CSV文件（除第一个文件外）并将其追加到合并后的文件中（跳过表头）
for file in "$folder_path"/*second_part.csv; do
    if [ -e "$file" ] ; then
        # 使用 tail 命令跳过CSV文件的标题行，并将其追加到合并后的文件中
        tail -n +2 "$file" >> "$output_file"
        # echo "文件 $file 已合并到 $output_file"
    fi
done

echo "合并完成。合并后的文件保存在: $output_file"
mv $output_file ../result



