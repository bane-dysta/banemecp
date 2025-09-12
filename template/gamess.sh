#!/bin/bash

# 环境变量定义
MOKIT_ENV='source ~/.bashrc
conda activate mokit'
G16_ENV='source /apps/gaussian16/env.sh'
GMS_ENV='export PATH=$PATH:/apps/gamess/gamess_src'

# 默认参数
MEM=96
NPROCS=36
FUNCTIONAL="BHANDHLYP"
BASIS="def2svp"
CHARGE=0
DFTTYP=""  # 空值表示自动映射
MWORDS=300
RUNTYP="GRADIENT"
IROOT=""
MULT=""   # 新增：多重度参数
CUSTOM_BASENAME=""

# 泛函映射表（Gaussian -> GAMESS）
declare -A FUNCTIONAL_MAP
FUNCTIONAL_MAP["M062X"]="M06-2X"
FUNCTIONAL_MAP["CAM-B3LYP"]="CAMB3LYP"
FUNCTIONAL_MAP["BHANDHLYP"]="BHHLYP"
FUNCTIONAL_MAP["BHandHLYP"]="BHHLYP"

# 函数：显示使用说明
show_usage() {
    echo "用法: $0 <xyz_file> [选项]"
    echo "选项："
    echo "  --mem <value>        内存大小 (默认: $MEM GB)"
    echo "  --nprocs <value>     CPU核心数 (默认: $NPROCS)"
    echo "  --functional <value> 泛函 (默认: $FUNCTIONAL)"
    echo "  --basis <value>      基组 (默认: $BASIS)"
    echo "  --dfttyp <value>     GAMESS泛函 (默认: 自动映射)"
    echo "  --mwords <value>     内存字数 (默认: $MWORDS)"
    echo "  --runtyp <value>     运行类型 ENERGY或GRADIENT (默认: $RUNTYP)"
    echo "  --iroot <value>      在TDDFT部分添加IROOT=n参数"
    echo "  --mult <value>       在TDDFT部分添加MULT=n参数（多重度）"
    echo "  --basename <value>   自定义输入输出文件名前缀 (默认: 使用xyz文件名)"
    echo "  -h, --help           显示此帮助信息"
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --mem)
            MEM="$2"
            shift 2
            ;;
        --nprocs)
            NPROCS="$2"
            shift 2
            ;;
        --functional)
            FUNCTIONAL="${2^^}"  # 转换为大写
            shift 2
            ;;
        --basis)
            BASIS="$2"
            shift 2
            ;;
        --dfttyp)
            DFTTYP="${2^^}"  # 转换为大写
            shift 2
            ;;
        --mwords)
            MWORDS="$2"
            shift 2
            ;;
        --runtyp)
            RUNTYP="${2^^}"  # 转换为大写
            shift 2
            ;;
        --iroot)
            IROOT="$2"
            shift 2
            ;;
        --mult)
            MULT="$2"
            shift 2
            ;;
        --basename)
            CUSTOM_BASENAME="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        -*)
            echo "未知选项: $1"
            show_usage
            exit 1
            ;;
        *)
            if [[ -z $XYZ_FILE ]]; then
                XYZ_FILE="$1"
            else
                echo "错误: 多余的参数 $1"
                show_usage
                exit 1
            fi
            shift
            ;;
    esac
done

# 检查是否提供了xyz文件
if [[ -z $XYZ_FILE ]]; then
    echo "错误: 请提供xyz文件"
    show_usage
    exit 1
fi

# 检查xyz文件是否存在
if [[ ! -f $XYZ_FILE ]]; then
    echo "错误: 文件 $XYZ_FILE 不存在"
    exit 1
fi

# 验证RUNTYP参数
if [[ "$RUNTYP" != "ENERGY" && "$RUNTYP" != "GRADIENT" ]]; then
    echo "错误: RUNTYP必须是ENERGY或GRADIENT"
    exit 1
fi

# 验证MULT参数（如果提供了的话）
if [[ -n $MULT ]] && ! [[ "$MULT" =~ ^[1-9][0-9]*$ ]]; then
    echo "错误: MULT必须是正整数"
    exit 1
fi

# 确定文件名前缀
if [[ -n $CUSTOM_BASENAME ]]; then
    BASENAME="$CUSTOM_BASENAME"
    echo "使用自定义文件名前缀: $BASENAME"
else
    BASENAME=$(basename "$XYZ_FILE" .xyz)
    echo "使用默认文件名前缀（基于xyz文件名）: $BASENAME"
fi

# 从xyz文件读取电荷
CHARGE=$(sed -n '2p' "$XYZ_FILE" | awk '{print $1}')
if ! [[ "$CHARGE" =~ ^-?[0-9]+$ ]]; then
    echo "警告: 无法从xyz文件读取电荷，使用默认值 0"
    CHARGE=0
fi

# 读取几何结构（跳过前两行）
GEOM=$(tail -n +3 "$XYZ_FILE")

# 处理泛函映射
if [[ -n ${FUNCTIONAL_MAP[$FUNCTIONAL]} ]]; then
    GAMESS_FUNCTIONAL=${FUNCTIONAL_MAP[$FUNCTIONAL]}
else
    GAMESS_FUNCTIONAL=$FUNCTIONAL
fi

# 如果用户没有手动指定DFTTYP，使用映射后的泛函
if [[ -z $DFTTYP ]]; then
    DFTTYP=$GAMESS_FUNCTIONAL
fi

# 打印当前使用的变量
echo "========================================="
echo "分子计算自动化脚本 - 参数列表"
echo "========================================="
printf "%-20s %-15s %-20s %-15s\n" "参数" "值" "参数" "值"
echo "-----------------------------------------"
printf "%-20s %-15s %-20s %-15s\n" "XYZ文件" "$XYZ_FILE" "内存(GB)" "$MEM"
printf "%-20s %-15s %-20s %-15s\n" "CPU核心数" "$NPROCS" "泛函" "$FUNCTIONAL"
printf "%-20s %-15s %-20s %-15s\n" "基组" "$BASIS" "电荷" "$CHARGE"
printf "%-20s %-15s %-20s %-15s\n" "GAMESS泛函" "$DFTTYP" "内存字数" "$MWORDS"
printf "%-20s %-15s %-20s %-15s\n" "运行类型" "$RUNTYP" "IROOT" "${IROOT:-未设置}"
printf "%-20s %-15s %-20s %-15s\n" "MULT" "${MULT:-未设置}" "文件名前缀" "$BASENAME"
echo "========================================="

# 步骤1: 生成G16输入文件
echo "步骤1: 生成G16输入文件..."
G16_INPUT="${BASENAME}.gjf"

cat > "$G16_INPUT" << EOF
%chk=${BASENAME}.chk
%mem=${MEM}GB
%nprocshared=${NPROCS}
#p RO${FUNCTIONAL}/${BASIS} nosymm int=nobasistransform

title

${CHARGE} 3
${GEOM}

EOF

echo "已生成: $G16_INPUT"

# 加载G16环境并运行计算
echo "步骤2: 运行G16计算..."
eval "$G16_ENV"
g16 "$G16_INPUT"

if [[ $? -ne 0 ]]; then
    echo "错误: G16计算失败"
    exit 1
fi

# 运行formchk
echo "步骤3: 运行formchk..."
formchk "${BASENAME}.chk" "${BASENAME}.fchk"

if [[ $? -ne 0 ]]; then
    echo "错误: formchk失败"
    exit 1
fi

# 步骤4: 加载MOKIT环境并处理fchk文件
echo "步骤4: 使用MOKIT处理fchk文件..."
eval "$MOKIT_ENV"
fch2inp "${BASENAME}.fchk" -mrsf

if [[ $? -ne 0 ]]; then
    echo "错误: fch2inp失败"
    exit 1
fi

# 步骤5: 修改生成的inp文件
echo "步骤5: 修改GAMESS输入文件..."
INP_FILE="${BASENAME}.inp"
if [[ ! -f $INP_FILE ]]; then
    echo "错误: 找不到生成的inp文件"
    exit 1
fi

# 备份原文件
cp "$INP_FILE" "${INP_FILE}.bak"

# 替换DFTTYP
sed -i "s/DFTTYP=[A-Z0-9-]*/DFTTYP=${DFTTYP}/g" "$INP_FILE"
echo "已将GAMESS输入文件中的泛函修改为: $DFTTYP"

# 替换RUNTYP
sed -i "s/RUNTYP=[A-Z]*/RUNTYP=${RUNTYP}/g" "$INP_FILE"
echo "已将RUNTYP修改为: $RUNTYP"

# 处理TDDFT参数（IROOT和MULT）
if [[ -n $IROOT || -n $MULT ]]; then
    # 构建要添加的TDDFT参数
    TDDFT_PARAMS=""
    if [[ -n $IROOT ]]; then
        TDDFT_PARAMS="${TDDFT_PARAMS} IROOT=${IROOT}"
    fi
    if [[ -n $MULT ]]; then
        TDDFT_PARAMS="${TDDFT_PARAMS} MULT=${MULT}"
    fi
    
    # 在$TDDFT行中添加参数
    sed -i "/\$TDDFT/s/NSTATE=[0-9]*/&${TDDFT_PARAMS}/" "$INP_FILE"
    
    if [[ -n $IROOT ]]; then
        echo "已在TDDFT部分添加IROOT=${IROOT}"
    fi
    if [[ -n $MULT ]]; then
        echo "已在TDDFT部分添加MULT=${MULT}"
    fi
fi

# 步骤6: 加载GMS环境并运行计算
echo "步骤6: 运行GAMESS计算..."
eval "$GMS_ENV"
rungms "$INP_FILE" 01 "$NPROCS" > "${BASENAME}.gms"

if [[ $? -ne 0 ]]; then
    echo "错误: GAMESS计算失败"
    exit 1
fi

echo "========================================="
echo "计算完成！"
echo "生成的文件："
echo "  - G16输入: $G16_INPUT"
echo "  - G16检查点: ${BASENAME}.chk"
echo "  - G16格式化检查点: ${BASENAME}.fchk"
echo "  - GAMESS输入: $INP_FILE"
echo "  - GAMESS输出: ${BASENAME}.gms"
echo "========================================"