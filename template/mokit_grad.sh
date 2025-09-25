#!/bin/bash
# ./mokit_grad.sh input_file [output_file]

input="$1"
output="${2:-$input}"  # 如果没有指定输出文件，就用输入文件

[ ! -f "$input" ] && { echo "文件不存在: $input"; exit 1; }

# 提取梯度并转换为三列格式，直接追加到文件
{
    echo ""
    echo "Gradients for baneMECP (gen by mokit_grad.sh)"
    
    # 提取所有梯度数值并重新格式化
    sed -n '/Cartesian gradients\|GRADIENT OF THE ENERGY/,/^$/p' "$input" | \
    grep -E '^[[:space:]]*[-+]?[0-9]' | \
    tr -s ' ' | tr ' ' '\n' | \
    grep -E '^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$' | \
    awk '{
        if(NR%3==1) printf "  %15.8E", $1
        else if(NR%3==2) printf "  %15.8E", $1  
        else printf "  %15.8E\n", $1
    }'
    echo ""
} >> "$output"



