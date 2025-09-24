#!/bin/bash

# 环境变量定义
MOKIT_ENV='source ~/.bashrc; conda activate mokit'
G16_ENV='source /apps/gaussian16/env.sh'
GMS_ENV='export PATH=$PATH:/apps/gamess/gamess_src'

# 默认参数
MEM=96; NPROCS=36; FUNCTIONAL="BHANDHLYP"; BASIS="def2svp"; CHARGE=0; DFTTYP=""
MWORDS=300; RUNTYP="GRADIENT"; IROOT=""; MULT=""; CUSTOM_BASENAME=""
FCHK_FILE=""  # RODFT波函数文件

# 泛函映射表（Gaussian -> GAMESS）
declare -A FUNCTIONAL_MAP
FUNCTIONAL_MAP["M062X"]="M06-2X"
FUNCTIONAL_MAP["CAM-B3LYP"]="CAMB3LYP"
FUNCTIONAL_MAP["BHANDHLYP"]="BHHLYP"
FUNCTIONAL_MAP["BHandHLYP"]="BHHLYP"

# 函数：显示使用说明
show_usage() {
    echo "用法: $0 <input_file> [选项]"
    echo "输入文件："
    echo "  xyz文件：完整流程（G16 → MOKIT → GAMESS）"
    echo "  fchk文件：从步骤4开始（MOKIT → GAMESS）"
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
    echo "  --basename <value>   自定义输入输出文件名前缀 (默认: 使用输入文件名)"
    echo "  -h, --help           显示此帮助信息"
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --mem) MEM="$2"; shift 2;;
        --nprocs) NPROCS="$2"; shift 2;;
        --functional) FUNCTIONAL="${2^^}"; shift 2;;
        --basis) BASIS="$2"; shift 2;;
        --dfttyp) DFTTYP="${2^^}"; shift 2;;
        --mwords) MWORDS="$2"; shift 2;;
        --runtyp) RUNTYP="${2^^}"; shift 2;;
        --iroot) IROOT="$2"; shift 2;;
        --mult) MULT="$2"; shift 2;;
        --basename) CUSTOM_BASENAME="$2"; shift 2;;
        -h|--help) show_usage; exit 0;;
        -*) echo "未知选项: $1"; show_usage; exit 1;;
        *) if [[ -z $INPUT_FILE ]]; then INPUT_FILE="$1"; else echo "错误: 多余的参数 $1"; show_usage; exit 1; fi; shift;;
    esac
done

# 检查输入文件和参数
[[ -z $INPUT_FILE ]] && { echo "错误: 请提供输入文件（xyz或fchk）"; show_usage; exit 1; }
[[ ! -f $INPUT_FILE ]] && { echo "错误: 文件 $INPUT_FILE 不存在"; exit 1; }
# 检查文件类型并设置相应变量
if [[ "$INPUT_FILE" == *.xyz ]]; then
    XYZ_FILE="$INPUT_FILE"
    START_STEP=1
    echo "检测到xyz文件，将执行完整流程（G16 → MOKIT → GAMESS）"
elif [[ "$INPUT_FILE" == *.fchk ]]; then
    FCHK_FILE="$INPUT_FILE"
    START_STEP=4
    echo "检测到fchk文件，将从步骤4开始执行（MOKIT → GAMESS）"
else
    echo "错误: 不支持的文件类型。请提供.xyz或.fchk文件"
    exit 1
fi

[[ "$RUNTYP" != "ENERGY" && "$RUNTYP" != "GRADIENT" ]] && { echo "错误: RUNTYP必须是ENERGY或GRADIENT"; exit 1; }
[[ -n $MULT ]] && ! [[ "$MULT" =~ ^[1-9][0-9]*$ ]] && { echo "错误: MULT必须是正整数"; exit 1; }

# 确定文件名前缀
if [[ -n $CUSTOM_BASENAME ]]; then
    BASENAME="$CUSTOM_BASENAME"
else
    if [[ -n $XYZ_FILE ]]; then
        BASENAME=$(basename "$XYZ_FILE" .xyz)
    else
        BASENAME=$(basename "$FCHK_FILE" .fchk)
    fi
fi

# 如果是xyz文件，读取电荷和几何结构
if [[ -n $XYZ_FILE ]]; then
    CHARGE=$(sed -n '2p' "$XYZ_FILE" | awk '{print $1}')
    [[ ! "$CHARGE" =~ ^-?[0-9]+$ ]] && { echo "警告: 无法从xyz文件读取电荷，使用默认值 0"; CHARGE=0; }
    GEOM=$(tail -n +3 "$XYZ_FILE")
fi

# 处理泛函映射
GAMESS_FUNCTIONAL=${FUNCTIONAL_MAP[$FUNCTIONAL]:-$FUNCTIONAL}
[[ -z $DFTTYP ]] && DFTTYP=$GAMESS_FUNCTIONAL

# 显示计算参数（简化版）
if [[ -n $XYZ_FILE ]]; then
    echo "计算参数: 文件=${XYZ_FILE} | 前缀=${BASENAME} | 泛函=${FUNCTIONAL}→${DFTTYP} | 基组=${BASIS} | 电荷=${CHARGE} | 内存=${MEM}GB | 核心=${NPROCS}"
else
    echo "计算参数: 文件=${FCHK_FILE} | 前缀=${BASENAME} | 泛函映射=${FUNCTIONAL}→${DFTTYP} | 内存=${MEM}GB | 核心=${NPROCS}"
fi
[[ -n $IROOT ]] && echo "TDDFT参数: IROOT=${IROOT}"
[[ -n $MULT ]] && echo "TDDFT参数: MULT=${MULT}"

# 步骤1-3: G16计算（仅当输入是xyz文件时执行）
if [[ $START_STEP -eq 1 ]]; then
    # 步骤1: 生成G16输入文件
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

    # 步骤2: 运行G16计算
    echo -n "运行G16计算..."
    eval "$G16_ENV"
    g16 "$G16_INPUT" > /dev/null 2>&1
    [[ $? -ne 0 ]] && { echo " 失败"; exit 1; } || echo " 完成"

    # 步骤3: 运行formchk
    echo -n "运行formchk..."
    formchk "${BASENAME}.chk" "${BASENAME}.fchk" > /dev/null 2>&1
    [[ $? -ne 0 ]] && { echo " 失败"; exit 1; } || echo " 完成"
    
    FCHK_FILE="${BASENAME}.fchk"
fi

# 步骤4: 使用MOKIT处理fchk文件
echo -n "MOKIT处理fchk文件..."
eval "$MOKIT_ENV"

# 如果basename与fchk文件名不同，需要先重命名fchk文件
ORIGINAL_FCHK="$FCHK_FILE"
if [[ "$(basename "$FCHK_FILE" .fchk)" != "$BASENAME" ]]; then
    NEW_FCHK_FILE="${BASENAME}.fchk"
    cp "$FCHK_FILE" "$NEW_FCHK_FILE"
    FCHK_FILE="$NEW_FCHK_FILE"
    echo " 重命名fchk文件: $(basename "$ORIGINAL_FCHK") → $(basename "$NEW_FCHK_FILE")"
fi

fch2inp "$FCHK_FILE" -mrsf > /dev/null 2>&1
[[ $? -ne 0 ]] && { echo " 失败"; exit 1; } || echo " 完成"

# 步骤5: 修改GAMESS输入文件
INP_FILE="${BASENAME}.inp"
[[ ! -f $INP_FILE ]] && { echo "错误: 找不到生成的inp文件"; exit 1; }

cp "$INP_FILE" "${INP_FILE}.bak"
sed -i "s/DFTTYP=[A-Z0-9-]*/DFTTYP=${DFTTYP}/g" "$INP_FILE"
sed -i "s/RUNTYP=[A-Z]*/RUNTYP=${RUNTYP}/g" "$INP_FILE"

# 处理TDDFT参数
if [[ -n $IROOT || -n $MULT ]]; then
    TDDFT_PARAMS=""
    [[ -n $IROOT ]] && TDDFT_PARAMS="${TDDFT_PARAMS} IROOT=${IROOT}"
    [[ -n $MULT ]] && TDDFT_PARAMS="${TDDFT_PARAMS} MULT=${MULT}"
    sed -i "/\$TDDFT/s/NSTATE=[0-9]*/&${TDDFT_PARAMS}/" "$INP_FILE"
fi

echo "GAMESS输入文件已修改: 泛函=${DFTTYP} | 运行类型=${RUNTYP}${TDDFT_PARAMS:+ | TDDFT参数=${TDDFT_PARAMS}}"

# 步骤6: 运行GAMESS计算
echo -n "运行GAMESS计算..."
eval "$GMS_ENV"
rungms "$INP_FILE" 01 "$NPROCS" > "${BASENAME}.gms"
if [[ $? -ne 0 ]]; then
    echo " 失败"
    exit 1
else
    echo " 完成"
fi

# 检查SCF是否收敛
if tail -100 "${BASENAME}.gms" | grep -q "SCF DID NOT CONVERGE"; then
    echo "检测到SCF未收敛，修改SCF参数后重新计算..."
    
    # 删除gms文件
    rm "${BASENAME}.gms"
    
    # 修改inp文件的SCF部分
    if grep -q "\$SCF.*DIRSCF=\.T\." "$INP_FILE"; then
        # 如果已经有DIRSCF=.T.，添加DIIS和SOSCF参数
        sed -i '/\$SCF.*DIRSCF=\.T\./s/\$END/ DIIS=.F. SOSCF=.T. $END/' "$INP_FILE"
    else
        # 如果没有找到DIRSCF=.T.，查找$SCF行并修改
        sed -i '/\$SCF/s/\$END/ DIRSCF=.T. DIIS=.F. SOSCF=.T. $END/' "$INP_FILE"
    fi
    
    echo "已修改SCF参数: DIRSCF=.T. DIIS=.F. SOSCF=.T."
    
    # 重新运行GAMESS计算
    echo -n "重新运行GAMESS计算..."
    rungms "$INP_FILE" 01 "$NPROCS" > "${BASENAME}.gms"
    if [[ $? -ne 0 ]]; then
        echo " 失败"
        exit 1
    else
        echo " 完成"
    fi
    
    # 再次检查是否收敛
    if tail -100 "${BASENAME}.gms" | grep -q "SCF DID NOT CONVERGE"; then
        echo "警告: SCF仍未收敛，请检查计算设置"
    else
        echo "SCF已收敛"
    fi
fi

echo "计算完成！生成文件: ${G16_INPUT:-无}, ${BASENAME}.chk${G16_INPUT:+, ${BASENAME}.fchk}, ${INP_FILE}, ${BASENAME}.gms"