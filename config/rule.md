# baneMECP 配置文件格式规则

配置文件采用INI格式，用于定义量子化学程序的执行环境和结果提取规则。

## 文件结构

```ini
[main]
# 主要配置参数

[state_type1]
# 第一种状态的提取规则

[state_type2]  
# 第二种状态的提取规则
```

## [main] 部分

### 必需参数

- `OUTPUT_SUFFIX`: 输出文件后缀名
  ```ini
  OUTPUT_SUFFIX=log  # 例如：step1_state1.log
  ```

### 可选参数

- `env`: 环境变量设置（多行字符串，用单引号包围）
  ```ini
  env='source ~/.bashrc
  module load gaussian
  export OMP_NUM_THREADS=8'
  ```

- `INPUT_SUFFIX`: 输入文件后缀名
  ```ini
  INPUT_SUFFIX=gjf  # 例如：step1_state1.gjf
  ```

- `RUN_CMD`: 程序执行命令（多行字符串，用单引号包围）
  ```ini
  RUN_CMD='echo "Starting calculation"
  g16 ${ACTUAL_INP} > ${ACTUAL_OUT}
  echo "Calculation done"'
  ```
  
  **占位符变量：**
  - `${ACTUAL_INP}`: 实际输入文件名
  - `${ACTUAL_OUT}`: 实际输出文件名

**注意：** 如果缺少`RUN_CMD`参数，程序将切换为external模式，直接执行`%InpTmplt`中的命令。

## 状态提取规则部分

每个状态类型对应一个section，section名称需要在`%GrpTmplt`中引用。

### 能量提取参数

- `E`: 能量提取的正则表达式
  ```ini
  E=SCF Done:\s+E\([^)]*\)\s*=\s*([-\d\.]+)
  ```
  - 正则表达式中的第一个捕获组`()`将被提取为能量值
  - 支持科学计数法和负数

- `E.LoacteCount`: 能量定位次数（可选，默认为1）
  ```ini
  E.LoacteCount=2  # 提取第2次出现的匹配
  ```

### 梯度提取参数

- `GARD.Loacte`: 梯度块的定位字符串
  ```ini
  GARD.Loacte=Forces (Hartrees/Bohr)
  ```

- `GARD.LoacteCount`: 梯度块定位次数（可选，默认为1）
  ```ini
  GARD.LoacteCount=2  # 定位到第2次出现的梯度块
  ```

- `GARD.NLineSkip`: 从定位字符串开始跳过的行数
  ```ini
  GARD.NLineSkip=2  # 跳过2行后开始读取梯度数据
  ```

- `GARD.TargetColumns`: 梯度数据的目标列（X,Y,Z坐标）
  ```ini
  GARD.TargetColumns=3,4,5  # 读取第3,4,5列作为X,Y,Z梯度
  ```

- `GARD.EndBy`: 梯度块的结束标记（可选）
  ```ini
  GARD.EndBy=----  # 遇到"----"时停止读取
  # 留空表示遇到空行时停止
  GARD.EndBy=
  ```

- `GARD.unit`: 梯度单位（可选）
  ```ini
  GARD.unit=hartree_per_bohr  # 目前仅支持hartree_per_bohr
  ```

- `GRAD.Type`: 梯度类型（可选，默认为force）
  ```ini
  GRAD.Type=grad   # 输出为梯度，需要取负号转换为力
  GRAD.Type=force  # 输出为力，无需转换
  ```

## 配置示例

### Gaussian 16配置
```ini
[main]
env='source ~/.bashrc
module load g16'
INPUT_SUFFIX=gjf
OUTPUT_SUFFIX=log
RUN_CMD='g16 ${ACTUAL_INP} > ${ACTUAL_OUT}'

[scf]
E=SCF Done:\s+E\([^)]*\)\s*=\s*([-\d\.]+)
GARD.Loacte=Forces (Hartrees/Bohr)
GARD.NLineSkip=2
GARD.TargetColumns=3,4,5
GARD.EndBy=----

[td]
E=Total Energy,\s*E\([^)]*\)\s*=\s*([-\d\.]+)
GARD.Loacte=Forces (Hartrees/Bohr)
GARD.NLineSkip=2
GARD.TargetColumns=3,4,5
GARD.EndBy=----
```

### GAMESS配置（external模式）
```ini
[main]
OUTPUT_SUFFIX=gms
# 注意：无RUN_CMD，将使用external模式

[mrsf1]
E=1\s{2}A\s+(-?\d+\.?\d+)
GARD.Loacte=GRADIENT OF THE ENERGY
GARD.LoacteCount=2
GARD.NLineSkip=3
GARD.TargetColumns=3,4,5
GARD.EndBy=
GRAD.Type=grad

[mrsf2]
E=2\s{2}A\s+(-?\d+\.?\d+)
GARD.Loacte=GRADIENT OF THE ENERGY
GARD.LoacteCount=2
GARD.NLineSkip=3
GARD.TargetColumns=3,4,5
GARD.EndBy=
GRAD.Type=grad
```

## 使用流程

1. **状态定义**: 在配置文件中定义各种电子状态的提取规则
2. **模板引用**: 在输入文件的`%GrpTmplt`中引用状态名称
3. **自动提取**: 程序根据规则自动从输出文件中提取能量和梯度
4. **单位转换**: 程序自动处理单位转换（Hartree/Bohr → Hartree/Angstrom）

