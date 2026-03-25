
```{=latex}
\begin{titlepage}
\centering
\vspace*{2.8cm}
{\Huge\bfseries baneMECP 软件说明书\par}
\vspace{0.9cm}
{\Large 面向双态交叉点搜索的迭代驱动、模板生成与结果提取系统\par}
\vspace{2.2cm}
\begin{tabular}{rl}
产品名称： & baneMECP \\
版本： & 1.1 \\
文档类型： & 软件说明书 \\
语言： & 中文 \\
\end{tabular}
\vfill
{\large Bane Project\par}
\vspace*{1.2cm}
\end{titlepage}

\tableofcontents
\newpage
```

# 概述

## 产品定位

baneMECP 是面向双态最低能量交叉点（MECP）搜索任务的迭代驱动程序。系统以标准 XYZ 几何、两态能量与梯度为核心输入，通过模板渲染、量子化学程序调用、结果提取和 MECP 求解器迭代更新，在统一的目录结构下完成从初始结构到最终交叉点结构的连续优化。

baneMECP 采用“驱动器 + 可配置接口 + 数值求解器”的工作方式组织任务：

- `banemecp`：负责读取输入文件、组织逐步迭代、渲染模板、调用外部程序、提取结果并管理工作目录。
- `<prog>.conf`：定义某一量子化学程序的运行命令、输入/输出后缀以及能量和梯度的提取规则。
- `baneMECP.f` / `MECP.x`：负责根据两态能量与梯度计算有效梯度、更新坐标并判断收敛。

## 适用场景

baneMECP 适用于以下场景：

- 两个电子态之间的最低能量交叉点搜索。
- 量子化学后端程序不固定，需要通过模板和配置文件接入 Gaussian、GAMESS、BDF、MOKIT 或站点自定义脚本。
- 两态计算需要在每一步根据当前几何重新生成输入文件。
- 输出文件的能量和梯度位置随程序不同而变化，需要通过正则和块定位规则做提取。
- 需要在本地工作目录或临时目录中保留逐步结构、状态文件、收敛信息和轨迹文件。

## 核心能力

- 使用统一输入文件描述程序类型、初始几何、控制参数、两态输入模板和结果提取映射。
- 支持 Harvey、BPUPD、Lagrange / Penalty Function 三类 MECP 有效梯度算法。
- 支持 BFGS、GDIIS、GEDIIS 三类几何更新加速策略。
- 支持 BFGS、PSB、SR1、Bofill 四类 Hessian 更新方式。
- 支持内置模式与 external 模式两类量子化学接口方式。
- 支持从量子化学输出中自动提取两态能量与梯度，也支持由外部脚本直接写出 `.grad1` / `.grad2`。
- 支持重启、临时目录运行、保留或清理中间文件、逐步轨迹记录和最终结构输出。

## 核心对象

| 对象 | 说明 |
| --- | --- |
| 输入文件 | 以 `%Prog`、`%Geom`、`%Control`、`%InpTmplt1`、`%InpTmplt2`、`%GrpTmplt` 组织的任务描述文件。 |
| 配置文件 | 形如 `g16.conf`、`gamess.conf` 的 INI 文件，用于定义程序运行方式和结果提取规则。 |
| 当前步 | 一次 MECP 迭代，对应 `astepN` 一组输入、输出、梯度和结构文件。 |
| `.grad1` / `.grad2` | 两个电子态的能量与梯度文件，是求解器的直接输入。 |
| `MECP.state` | 保存上一步几何、梯度、逆 Hessian 和算法状态，用于重启。 |
| `convg.tmp` | 保存当前步的收敛状态、能量差、梯度与位移指标。 |
| `opt.trj` | 逐步结构轨迹文件，第二行写入收敛摘要。 |
| `<geom_stem>_final.xyz` | 收敛后写出的最终结构。 |

# 系统组成与运行模型

## 组成部分

当前版本由以下三个层次组成：

### 迭代驱动层

`iter/banemecp.cpp` 编译生成命令行程序 `banemecp`。该程序负责：

- 解析输入文件；
- 读取 `%Control` 控制参数；
- 解析 `%GrpTmplt` 指向的状态类型；
- 搜索并读取 `<prog>.conf`；
- 搜索并复制 `baneMECP.f`；
- 按原子数修改 `baneMECP.f` 并编译 `MECP.x`；
- 为每一步生成两态输入；
- 顺序执行 state1 与 state2 计算；
- 根据配置提取能量和梯度，或直接接收外部脚本生成的 `.grad` 文件；
- 调用 `MECP.x` 计算新几何；
- 判断收敛、写出最终结果并执行清理。

### 程序接口层

接口层由 `%InpTmplt1` / `%InpTmplt2` 与 `<prog>.conf` 共同构成：

- `%InpTmplt1` / `%InpTmplt2` 定义逐步输入模板或逐步执行脚本模板；
- `<prog>.conf` 的 `[main]` 段定义环境与执行命令；
- `<prog>.conf` 的状态段定义当前程序的能量和梯度提取方式。

### 求解器层

`opt/baneMECP.f` 在运行时被复制并改写为 `MECP_temp.f`，再使用 `gfortran` 编译为 `MECP.x`。求解器负责：

- 读取 `astepN.xyz`、`astepN.grad1`、`astepN.grad2`；
- 计算有效梯度；
- 更新几何和逆 Hessian；
- 根据所选算法判断是否收敛；
- 写出 `new.xyz` 或 `final.xyz`；
- 维护 `MECP.state`、`convg.tmp` 和 `MECP.diis`。

## 标准处理流程

baneMECP 的标准运行流程如下：

1. 读取输入文件并解析 `%Prog`、`%Geom`、`%Control`、模板段和 `%GrpTmplt`。
2. 解析控制参数，确定 `maxcyc`、`tmpdir`、`restart`、算法、阈值及加速参数。
3. 以当前工作目录为基准搜索 `<prog>.conf` 与 `baneMECP.f`。
4. 如设置 `tmpdir` 且其值不为 `.`，创建并切换到临时目录。
5. 复制几何文件、`baneMECP.f` 以及选中的配置文件到工作目录。
6. 根据初始 XYZ 原子数改写 `baneMECP.f` 中的 `Natom` 常量并编译 `MECP.x`。
7. 生成 `astep1.xyz`，并将其写入轨迹文件 `opt.trj`。
8. 对每个迭代步 `N`：
   1. 按模板渲染 `astepN_state1.*` 与 `astepN_state2.*`；
   2. 顺序执行两态计算；
   3. 生成或读取 `astepN.grad1` 与 `astepN.grad2`；
   4. 调用 `./MECP.x astepN ...`；
   5. 读取 `convg.tmp` 判断是否收敛；
   6. 若未收敛，将 `new.xyz` 重命名为 `astep(N+1).xyz` 并进入下一步；
   7. 若收敛，写出 `<geom_stem>_final.xyz` 并结束。
9. 根据 `keeptmp` 与当前目录模式决定保留或清理中间文件。

## 相对路径与搜索顺序

### 当前工作目录的作用

当前版本的相对路径解析基准是**程序启动时的当前工作目录**，而不是输入文件所在目录。以下对象均按当前工作目录搜索或解析：

- `%Geom` 指定的几何文件；
- `<prog>.conf`；
- `baneMECP.f`；
- `tmpdir` 相对路径；
- 模板脚本中自行引用的本地文件。

因此，推荐在任务目录中直接执行 `banemecp input.inp`。

### `baneMECP.f` 搜索顺序

`baneMECP.f` 按以下顺序搜索：

1. 当前工作目录
2. `banemecp` 可执行文件所在目录
3. `~/.bane/mecp/`

### `<prog>.conf` 搜索顺序

对于 `%Prog g16` 之类的程序名，配置文件 `g16.conf` 按以下顺序搜索：

1. 当前工作目录
2. `banemecp` 可执行文件所在目录
3. `~/.bane/mecp/`

## 推荐目录布局

推荐将程序、求解器源文件和配置文件组织为如下布局：

```text
~/.bane/mecp/
  banemecp
  baneMECP.f
  g16.conf
  gamess.conf
  bdf.conf
  mokit.conf
```

单个任务目录建议如下：

```text
case01/
  run.inp
  mol.xyz
  gamess.sh
  tmp/
```

当模板依赖额外脚本或辅助文件时，这些文件应放置在当前工作目录，或在模板脚本中以显式路径复制到 `tmpdir`。

# 安装与构建

## 运行环境

当前版本面向类 Unix 环境设计，运行时使用以下系统能力：

- `bash`：external 模式固定以 `bash` 执行逐步脚本；
- `gfortran`：用于在运行时编译 `MECP.x`；
- 标准文件系统操作：创建目录、复制文件、删除文件、重命名文件；
- 已安装并可调用的量子化学程序或站点脚本。

## 构建 `banemecp`

项目提供的编译脚本位于 `iter/compile.sh`。典型构建方式如下：

```bash
cd iter
bash compile.sh
```

脚本当前包含两条编译命令：

```bash
g++ -std=c++17 -O2 banemecp.cpp -o banemecp -lstdc++fs
g++ -std=c++17 -O3 -DNDEBUG -static-libgcc -static-libstdc++ -s banemecp.cpp -o banemecp -lstdc++fs
```

构建完成后，将生成可执行文件 `banemecp`。

## 运行时依赖

在执行 `banemecp` 前，应保证以下依赖已可用：

| 依赖 | 作用 |
| --- | --- |
| `gfortran` | 将 `baneMECP.f` 编译为 `MECP.x`。 |
| `%Prog` 对应的量子化学程序 | 执行 state1 / state2 两态计算。 |
| 对应程序的配置文件 | 内置模式与内置提取模式都依赖 `<prog>.conf`。 |
| 模板中调用的辅助脚本 | 如 `gamess.sh`、`mokit_grad.sh` 等站点脚本。 |

## 随包示例资源

源码包当前附带如下示例资源：

| 路径 | 作用 |
| --- | --- |
| `config/g16.conf` | Gaussian 16 示例配置。 |
| `config/gamess.conf` | GAMESS / MRSF 外部执行示例配置。 |
| `config/bdf.conf` | BDF 示例配置。 |
| `config/mokit.conf` | MOKIT 示例配置。 |
| `template/g16.inp` | Gaussian 输入文件示例。 |
| `template/gamess.inp` | GAMESS 外部脚本示例。 |
| `template/bdf.inp` | BDF 输入文件示例。 |
| `template/mokit.inp` | MOKIT 示例。 |
| `template/gamess.sh` | GAMESS 辅助脚本示例。 |
| `template/mokit_grad.sh` | MOKIT 梯度整理脚本示例。 |

这些文件用于站点接入参考。`env`、`RUN_CMD` 和脚本中的路径应按本地环境配置。

# 命令行用法

## 主命令

`banemecp` 当前只有一个用户入口形式：

```text
banemecp <inputfile.inp>
```

示例：

```bash
banemecp g16.inp
banemecp run.inp
```

## 命令行参数范围

当前版本不提供 `--help`、`--version`、`--debug` 之类的驱动级命令行选项。所有运行控制均通过输入文件中的 `%Control` 段指定。

## 启动后的内部动作

命令启动后，程序将自动完成以下内部动作：

- 解析输入文件；
- 编译与当前原子数匹配的 `MECP.x`；
- 逐步生成量子化学输入或脚本；
- 顺序执行两态计算；
- 组织 `convg.tmp`、`MECP.state`、`opt.trj` 与最终结构输出。

# 输入文件与 DSL 语法

## 文件结构

输入文件由若干指令段组成，段名以 `%` 开头。当前正式支持以下指令：

| 指令 | 作用 |
| --- | --- |
| `%Prog <prog>` | 指定程序名称并决定配置文件搜索名 `<prog>.conf`。 |
| `%Geom <file.xyz>` | 指定初始几何文件。 |
| `%Control` | 指定最大步数、临时目录、算法、阈值和加速参数。 |
| `%InpTmplt1` | 第 1 态输入模板或脚本模板。 |
| `%InpTmplt2` | 第 2 态输入模板或脚本模板。 |
| `%GrpTmplt` | 指定 state1 与 state2 的结果提取规则名称。 |

除模板正文外，输入文件不定义独立的注释语法。控制段中的值建议单独成行，不与注释混写。

## 词法与大小写规则

- 指令名大小写不敏感，例如 `%Prog` 与 `%prog` 等价。
- `%Control` 中的键名按小写处理，`algorithm=Harvey`、`algorithm=harvey` 等价。
- `%Prog` 的程序名按原样使用，文件名搜索区分大小写。推荐统一使用小写，如 `g16`、`gamess`、`bdf`、`mokit`、`external`。
- 段结束符必须写为单独一行 `end`。当前实现按修剪空白后的文本与 `end` 精确比较，`END` 不作为结束标志。

## `%Prog`

```text
%Prog g16
%Prog gamess
%Prog bdf
%Prog mokit
%Prog external
```

### 作用

`%Prog` 同时决定以下行为：

- 配置文件搜索名 `<prog>.conf`；
- 内置模式或 external 模式的判定起点；
- `%GrpTmplt` 引用的状态段所属配置文件。

### 模式判定规则

baneMECP 采用以下规则确定运行模式：

| 条件 | 结果 |
| --- | --- |
| `%Prog external` | 直接进入 external 模式。 |
| `%Prog <prog>` 且找到 `<prog>.conf`， | 内置模式。 |
| 并且 `[main]` 段含 `RUN_CMD` |  |
| `%Prog <prog>` 且找到 `<prog>.conf`， | external 模式。 |
| 但 `[main]` 段不含 `RUN_CMD` | |
| `%Prog <prog>` 且找不到 `<prog>.conf` | 报错退出。`external` 除外。 |

### 推荐理解方式

- **内置模式**：模板生成输入文件，程序根据 `RUN_CMD` 调用量子化学程序。
- **external 模式**：模板直接渲染为脚本并以 `bash` 执行；脚本可以自行完成量子化学计算、后处理和 `.grad` 文件生成。

对于需要保留程序名并使用内置提取规则的 external 工作流，推荐使用实际程序名（如 `gamess`）并在 `<prog>.conf` 中省略 `RUN_CMD`，而不是使用 `%Prog external`。

## `%Geom`

```text
%Geom mol.xyz
```

### 作用

`%Geom` 指定初始几何。当前版本要求输入为标准 XYZ 文件。

### 文件要求

XYZ 文件应满足以下格式：

1. 第 1 行：原子数；
2. 第 2 行：注释行；
3. 后续每行：`Element  X  Y  Z`。

示例：

```text
3
CH2 initial geometry
C   0.000000   0.000000   0.000000
H   0.000000   0.000000   1.080000
H   1.000000   0.000000  -0.360000
```

### 结果命名

最终收敛结构文件名由 `%Geom` 所指定文件的 stem 决定：

- `%Geom ch2.xyz` -> 输出 `ch2_final.xyz`
- `%Geom path/to/FeO.xyz` -> 输出 `FeO_final.xyz`

## `%Control`

`%Control` 段由若干 `key=value` 组成，以 `end` 结束。当前版本仅解析下表所列键。

### 基本控制参数

| 键 | 类型 | 默认值 | 说明 |
| --- | --- | --- | --- |
| `maxcyc` | 整数 | `50` | 最大迭代步数。 |
| `tmpdir` | 字符串 | `.` | 工作目录；`.` 表示在当前目录直接运行。 |
| `keeptmp` | 布尔 | `false` | 是否保留临时目录或当前目录中的中间文件。 |
| `debug` | 布尔 | `false` | 是否输出调试信息。 |
| `restart` | 布尔 | `false` | 是否从已有 `MECP.state` / `convg.tmp` 继续。 |

### 标准收敛阈值

| 键 | 默认值 | 单位 | 说明 |
| --- | --- | --- | --- |
| `tde` | `5e-5` | Hartree | 两态能量差阈值。 |
| `tgmax` | `7e-4` | Hartree/Angstrom | 最大梯度阈值。 |
| `tgrms` | `5e-4` | Hartree/Angstrom | RMS 梯度阈值。 |
| `tdxmax` | `4e-3` | Angstrom | 最大位移阈值。 |
| `tdxrms` | `2.5e-3` | Angstrom | RMS 位移阈值。 |
| `stpmx` | `0.1` | Angstrom | 信赖半径 / 最大步长。 |

### MECP 算法与优化方法

| 键 | 默认值 | 可取值 | 说明 |
| --- | --- | --- | --- |
| `algorithm` | `harvey` | `harvey`、`bpupd`、`lagrange`、`pf`、`lc` | 选择有效梯度算法。`lc` 在解析时映射为 `lagrange`。 |
| `opt_method` | `bfgs` | `bfgs`、`gdiis`、`gediis` | 几何更新加速方法。 |
| `ndiis` | `4` | 正整数 | DIIS 历史长度。求解器内部上限为 `8`。 |
| `hupd_method` | `bfgs` | `bfgs`、`psb`、`sr1`、`bofill` | Hessian 更新方式。 |

### Lagrange / PF 参数

以下参数在 `algorithm=lagrange` 或 `algorithm=pf` 时生效：

| 键 | 默认值 | 说明 |
| --- | --- | --- |
| `pf_alpha` | `0.02` | 罚函数参数 `alpha`。 |
| `pf_sigma` | `3.50` | 罚函数参数 `sigma`。 |
| `pf_tstep` | `1e-6` | PF 目标函数变化阈值。 |
| `pf_tgrad` | `5e-3` | PF 梯度阈值。 |
| `pf_thresh` | 无 | 同时设置 `pf_tstep` 与 `pf_tgrad`，接受 `1e-6,0.005` 或 `[1e-6, 0.005]`。 |
| `pf_thresh_step` | 无 | 等价于单独设置 `pf_tstep`。 |
| `pf_thresh_grad` | 无 | 等价于单独设置 `pf_tgrad`。 |

### 兼容键别名

`%Control` 当前还接受以下别名：

| 规范键 | 可接受别名 |
| --- | --- |
| `opt_method` | `opt`、`optmethod`、`diis_method` |
| `ndiis` | `diis_size`、`diissize` |
| `hupd_method` | `hupd`、`hess_update`、`hessupd` |

另外，`diis=<value>` 也会影响 `opt_method`：

- `diis=none`、`off`、`false` -> 采用 `bfgs`
- `diis=diis` -> 采用 `gdiis`
- `diis=gdiis` 或 `gediis` -> 采用对应方法
- 其他非空值 -> 采用 `gdiis`

### 布尔值规则

布尔键仅将字符串 `true` 识别为真；当前实现中，其他值都不会被当作真值处理。推荐写法如下：

```text
keeptmp=true
debug=false
restart=false
```

## `%InpTmplt1` 与 `%InpTmplt2`

### 作用

- `%InpTmplt1`：定义 state1 的输入模板或脚本模板。
- `%InpTmplt2`：定义 state2 的输入模板或脚本模板。

在标准两态 MECP 工作流中，两个模板都应提供。当前程序在以下场景中允许缺省 `%InpTmplt2`：

- external 模式；
- 且不使用 `%GrpTmplt` 内置提取；
- 且 `%InpTmplt1` 渲染后的脚本能自行生成当前步所需的 `astepN.grad1` 和 `astepN.grad2`。

### 内置模式下的生成物

在内置模式下，模板渲染后生成：

- `astepN_state1.<INPUT_SUFFIX>`
- `astepN_state2.<INPUT_SUFFIX>`

其中 `INPUT_SUFFIX` 由 `[main]` 段的 `INPUT_SUFFIX` 指定，例如 `gjf`、`inp`。

### external 模式下的生成物

在 external 模式下，模板渲染后生成：

- `astepN_state1.sh`
- `astepN_state2.sh`

程序随后执行：

```bash
bash astepN_state1.sh
bash astepN_state2.sh
```

### 模板占位符

当前版本在模板渲染时支持以下占位符：

| 占位符 | 含义 |
| --- | --- |
| `[name]` | 当前步基名，如 `astep1`、`astep2`。 |
| `[xyzfile]` | 当前步坐标文件名，如 `astep1.xyz`。 |
| `[geom]` | 当前步几何正文，格式为 `Element x y z`。 |
| `[nstep]` | 当前步号。 |
| `[tmpdir]` | `%Control` 中 `tmpdir` 的原始字符串值。 |

### 条件方括号规则

模板中的剩余方括号片段会按以下规则处理：

- 第 1 步：任意仍保留为 `[text]` 形式的片段会整体删除；
- 第 2 步及以后：`[text]` 会替换为 `text`。

该规则通常用于从第二步开始启用的关键词，例如：

```text
# B3LYP/def2SVP force [guess(read)]
```

渲染结果为：

- 第 1 步：`# B3LYP/def2SVP force `
- 第 2 步及以后：`# B3LYP/def2SVP force guess(read)`

### 模板执行顺序

每一步中，程序始终先执行 state1，再执行 state2。两态计算之间没有内置并行调度。

## `%GrpTmplt`

### 作用

`%GrpTmplt` 用于将 state1 / state2 映射到配置文件中的状态段，启用内置结果提取。

标准写法如下：

```text
%GrpTmplt
state1=scf
state2=td
end
```

### 解析规则

当前实现对 `%GrpTmplt` 的解析规则如下：

- 只读取段中的前两行；
- 第 1 行匹配 `state1=<type>`；
- 第 2 行匹配 `state2=<type>`；
- `<type>` 需能被正则 `\w+` 匹配，即通常由字母、数字或下划线组成；
- `<type>` 必须与 `<prog>.conf` 中的 section 名完全一致。

### 使用方式

- **提供 `%GrpTmplt`**：程序从输出文件中提取两态能量和梯度，并生成 `astepN.grad1` / `astepN.grad2`。
- **省略 `%GrpTmplt`**：仅适用于 external 模式，且模板脚本必须自行生成 `.grad1` / `.grad2`。

### 当前版本的实际要求

当使用内置提取时，程序总是尝试提取 state1 与 state2。内置提取工作流因此要求：

- 两态输出文件都已生成；
- `%GrpTmplt` 中同时给出 `state1` 和 `state2`；
- `%InpTmplt1` 和 `%InpTmplt2` 对应的两态计算均已完成。

# 配置文件格式

## 文件类型与结构

配置文件采用 INI 结构，通常命名为 `<prog>.conf`。当前版本识别以下部分：

1. `[main]`：运行命令与文件后缀；
2. 若干状态段：定义能量与梯度提取规则。

示例结构如下：

```ini
[main]
env='source ~/.bashrc
module load g16'
INPUT_SUFFIX=gjf
OUTPUT_SUFFIX=log
RUN_CMD='g16 ${ACTUAL_INP} > ${ACTUAL_OUT}'

[scf]
E=SCF Done:\s+E\([^)]*\)\s*=\s*([-\d\.]+)
GRAD.Loacte=Forces (Hartrees/Bohr)
GRAD.NLineSkip=2
GRAD.TargetColumns=3,4,5
GRAD.EndBy=----
GRAD.Type=force

[td]
E=Total Energy,\s*E\([^)]*\)\s*=\s*([-\d\.]+)
GRAD.Loacte=Forces (Hartrees/Bohr)
GRAD.NLineSkip=2
GRAD.TargetColumns=3,4,5
GRAD.EndBy=----
GRAD.Type=force
```

## 注释与多行值

### 注释

配置文件支持以 `#` 开头的整行注释。

### 多行值

`[main]` 段中的 `env` 和 `RUN_CMD` 可使用单引号括起的多行字符串：

```ini
env='source ~/.bashrc
module load g16'
```

解析规则如下：

- 值以首个单引号 `'` 开始；
- 持续读取到出现闭合单引号的行为止；
- 读取过程中保留换行。

## `[main]` 段

当前版本直接读取以下键：

| 键 | 是否必需 | 说明 |
| --- | --- | --- |
| `env` | 否 | 在执行 `RUN_CMD` 前先执行的环境初始化内容。 |
| `INPUT_SUFFIX` | 内置模式必需 | 输入文件后缀。 |
| `OUTPUT_SUFFIX` | 内置提取模式必需 | 输出文件后缀。 |
| `RUN_CMD` | 否 | 程序运行命令模板；若缺失，则驱动器切换为 external 模式。 |

### `env`

`env` 作为多行 shell 片段执行，常用于：

- `source ~/.bashrc`
- `module load g16`
- `export VAR=value`

### `RUN_CMD`

`RUN_CMD` 当前只支持两个占位符：

| 占位符 | 含义 |
| --- | --- |
| `${ACTUAL_INP}` | 当前步当前态输入文件名。 |
| `${ACTUAL_OUT}` | 当前步当前态输出文件名。 |

示例：

```ini
RUN_CMD='echo "Start g16"
g16 ${ACTUAL_INP} > ${ACTUAL_OUT}
echo "Done"'
```

### `RUN_CMD` 缺失时的行为

当配置文件存在但 `[main]` 段未提供 `RUN_CMD` 时，程序仍加载该配置文件，但实际运行模式切换为 external。该方式适用于：

- 使用脚本模板执行量子化学程序；
- 同时继续使用该配置文件中的状态段做内置结果提取。

## 状态段

状态段名称由 `%GrpTmplt` 引用。每个状态段描述一类输出文件中的能量和梯度提取规则。

### 当前版本实际读取的键

当前驱动程序对状态段读取以下键名：

| 键 | 是否必需 | 说明 |
| --- | --- | --- |
| `E` | 是 | 能量正则表达式；使用第一个捕获组作为能量值。 |
| `E.LoacteCount` | 否 | 选择第几次匹配的能量。默认 `1`。 |
| `GRAD.Loacte` | 是 | 梯度块定位字符串。 |
| `GRAD.LoacteCount` | 否 | 选择第几次出现的梯度定位字符串。默认 `1`。 |
| `GRAD.NLineSkip` | 是 | 从定位行之后额外跳过的行数。 |
| `GRAD.TargetColumns` | 是 | 梯度数据所在列，1 基计数，通常为三列。 |
| `GRAD.EndBy` | 是 | 梯度块提前结束判据；可为空字符串。 |
| `GRAD.Type` | 否 | `grad` 或 `force`，默认 `force`。 |

### `E`

`E` 应为一个能够提取能量数值的正则表达式。程序使用**第一个捕获组**作为能量字符串，并将其转换为浮点数。

示例：

```ini
E=SCF Done:\s+E\([^)]*\)\s*=\s*([-\d\.]+)
```

### `E.LoacteCount`

当同一输出中有多次能量匹配时，`E.LoacteCount` 指定使用第几次匹配结果：

```ini
E.LoacteCount=2
```

### `GRAD.Loacte`

`GRAD.Loacte` 用于定位梯度块开始处附近的文本。例如：

```ini
GRAD.Loacte=Forces (Hartrees/Bohr)
```

### `GRAD.NLineSkip`

从定位字符串所在行起，再跳过 `GRAD.NLineSkip` 行后，开始按原子数逐行读取梯度值。

### `GRAD.TargetColumns`

`GRAD.TargetColumns` 为 1 基列号，格式为逗号分隔的三列，例如：

```ini
GRAD.TargetColumns=3,4,5
```

表示从每一行的第 3、4、5 列读取 `x`、`y`、`z` 分量。

### `GRAD.EndBy`

`GRAD.EndBy` 用于检测梯度块是否过早结束：

- 非空时，若读取梯度行时遇到包含该字符串的行，则视为梯度块提前结束并报错；
- 为空时，若读取梯度行时遇到空白行，则视为梯度块提前结束并报错。

当前实现按原子数固定读取梯度行数，`GRAD.EndBy` 用于错误检查，而不是可变长度块提取。

### `GRAD.Type`

`GRAD.Type` 控制符号处理：

| 取值 | 行为 |
| --- | --- |
| `grad` | 按原样读取数值，不改符号。 |
| `force` | 对读取值整体乘以 `-1`，以转换为梯度方向。 |

该规则适用于输出“力”而不是“梯度”的程序。

### `GRAD.unit`

当前分发配置中常见 `GRAD.unit` 字段，但驱动器并不依据该键切换单位换算路径。当前结果链路假定输入给求解器的 `.grad` 数值来源为 Hartree/Bohr，求解器随后统一换算为 Hartree/Angstrom。

## 配置示例

### Gaussian 16

```ini
[main]
env='source ~/.bashrc
module load g16'
INPUT_SUFFIX=gjf
OUTPUT_SUFFIX=log
RUN_CMD='g16 ${ACTUAL_INP} > ${ACTUAL_OUT}'

[scf]
E=SCF Done:\s+E\([^)]*\)\s*=\s*([-\d\.]+)
GRAD.Loacte=Forces (Hartrees/Bohr)
GRAD.NLineSkip=2
GRAD.TargetColumns=3,4,5
GRAD.EndBy=----
GRAD.Type=force

[td]
E=Total Energy,\s*E\([^)]*\)\s*=\s*([-\d\.]+)
GRAD.Loacte=Forces (Hartrees/Bohr)
GRAD.NLineSkip=2
GRAD.TargetColumns=3,4,5
GRAD.EndBy=----
GRAD.Type=force
```

### external + built-in 提取

```ini
[main]
OUTPUT_SUFFIX=gms

[mrsf1]
E=1\s{2}A\s+(-?\d+\.?\d+)
GRAD.Loacte=GRADIENT OF THE ENERGY
GRAD.LoacteCount=2
GRAD.NLineSkip=3
GRAD.TargetColumns=3,4,5
GRAD.EndBy=
GRAD.Type=grad

[mrsf2]
E=2\s{2}A\s+(-?\d+\.?\d+)
GRAD.Loacte=GRADIENT OF THE ENERGY
GRAD.LoacteCount=2
GRAD.NLineSkip=3
GRAD.TargetColumns=3,4,5
GRAD.EndBy=
GRAD.Type=grad
```

上述配置不提供 `RUN_CMD`，因此驱动器进入 external 模式，但仍可借助 `%GrpTmplt` 做内置提取。

# 程序模式与输入生成

## 内置模式

### 判定条件

满足以下条件时进入内置模式：

- `%Prog <prog>`；
- 找到 `<prog>.conf`；
- `[main]` 段包含 `RUN_CMD`。

### 生成与执行流程

内置模式下，每一步按如下方式工作：

1. 渲染 `astepN_state1.<INPUT_SUFFIX>`；
2. 渲染 `astepN_state2.<INPUT_SUFFIX>`；
3. 依次执行 state1 与 state2 的 `RUN_CMD`；
4. 从 `astepN_state1.<OUTPUT_SUFFIX>` 和 `astepN_state2.<OUTPUT_SUFFIX>` 中提取结果；
5. 生成 `astepN.grad1` 与 `astepN.grad2`。

### 适用情况

- Gaussian、BDF 等能够直接通过统一命令执行的程序；
- 站点环境中已经有稳定的模块或环境初始化方式；
- 输出文件中能量和梯度位置稳定，可由配置文件描述。

## external 模式

### 判定条件

满足以下任一条件即为 external 模式：

- `%Prog external`；
- 或 `%Prog <prog>` 找到 `<prog>.conf`，但 `[main]` 段不含 `RUN_CMD`。

### 生成与执行流程

external 模式下，`%InpTmplt1` / `%InpTmplt2` 渲染为 `.sh` 脚本，并由程序使用 `bash` 直接执行。该模式适用于：

- 一步脚本中需要完成多程序串联；
- 需要在脚本内自定义文件复制、格式转换、后处理；
- 量子化学程序本身没有简单统一的单命令调用形式。

### 两种 external 工作方式

external 模式支持两类结果组织方式：

1. **external + 内置提取**
   - 脚本生成标准输出文件，如 `astepN_state1.gms`、`astepN_state2.gms`；
   - `%GrpTmplt` 提供 `state1` / `state2` 映射；
   - 配置文件中的状态段负责提取能量和梯度。

2. **external + 脚本直接产出 `.grad`**
   - 省略 `%GrpTmplt`；
   - 脚本自行写出 `astepN.grad1` 与 `astepN.grad2`；
   - 驱动器不再解析程序输出。

## 当前步输入命名规则

在所有模式下，步名统一为：

```text
astep1
astep2
astep3
...
```

相应文件命名规则如下：

| 文件 | 命名 |
| --- | --- |
| 当前步结构 | `astepN.xyz` |
| state1 输入 | `astepN_state1.<suffix>` 或 `astepN_state1.sh` |
| state2 输入 | `astepN_state2.<suffix>` 或 `astepN_state2.sh` |
| state1 梯度 | `astepN.grad1` |
| state2 梯度 | `astepN.grad2` |

# 算法、加速与收敛判据

## `algorithm`

### `harvey`

`harvey` 为默认算法，使用标准有效梯度方法，并采用以下五个判据联合判断收敛：

- 能量差；
- 最大梯度；
- RMS 梯度；
- 最大位移；
- RMS 位移。

### `bpupd`

`bpupd` 使用 BPUPD 有效梯度形式，并在 `MECP.state` 中额外保存 BPUPD 所需的历史向量。

### `lagrange` / `pf`

`lagrange` 与 `pf` 均走罚函数 / 增广拉格朗日路径。驱动器会向求解器传入：

- `pf_alpha`
- `pf_sigma`
- `pf_tstep`
- `pf_tgrad`

收敛时不再以标准五阈值为主，而是以 PF 三个函数判据为准。

## `opt_method`

### `bfgs`

- 默认方法；
- 使用 BFGS 类几何更新；
- 不启用 DIIS 历史文件加速。

### `gdiis`

- 在几何更新后追加 GDIIS 加速；
- DIIS 仅从第 3 步开始生效；
- 使用 `MECP.diis` 保存历史。

### `gediis`

- 在 GDIIS 基础上进一步尝试 GEDIIS；
- 若 GEDIIS 求解失败，则自动回退为 GDIIS 结果；
- 对 GDIIS 与 GEDIIS 提案做 0.5 权重混合。

## `ndiis`

`ndiis` 决定 DIIS 历史长度。当前求解器内部上限固定为 `8`：

- 小于 `1` 时自动钳制到 `1`；
- 大于 `8` 时自动钳制到 `8` 并输出警告。

## `hupd_method`

当前可选 Hessian 更新方式如下：

| 取值 | 说明 |
| --- | --- |
| `bfgs` | 默认 Hessian 更新。 |
| `psb` | Powell-symmetric-Broyden。 |
| `sr1` | Symmetric Rank-1。 |
| `bofill` | Bofill 更新。 |

## 标准收敛判据

在 `harvey` 与 `bpupd` 下，求解器使用以下五项标准联合判定：

| 指标 | 控制键 | 含义 |
| --- | --- | --- |
| `Energy_Gap` | `tde` | 两态能量差的绝对值。 |
| `Max_Gradient` | `tgmax` | 有效梯度最大分量。 |
| `RMS_Gradient` | `tgrms` | 有效梯度均方根。 |
| `Max_Displacement` | `tdxmax` | 当前步最大位移。 |
| `RMS_Displacement` | `tdxrms` | 当前步位移均方根。 |

五项全部满足时，`convg.tmp` 第一行写入 `CONVERGED`。

## PF / Lagrange 收敛判据

在 `lagrange` / `pf` 下，求解器使用以下三项判据：

| 指标 | 控制键 | 含义 |
| --- | --- | --- |
| `PF_Function1` | `pf_tstep` | 目标函数变化量。 |
| `PF_Function2` | `pf_tgrad` | PF 梯度函数 2。 |
| `PF_Function3` | `pf_tgrad` | PF 梯度函数 3。 |

三项全部满足时，视为收敛。

## `stpmx`

`stpmx` 控制每一步的信赖半径 / 最大更新步长。其作用范围包括：

- 标准几何更新；
- DIIS 步长控制。

较小的 `stpmx` 会减小几何移动幅度。

## `%Control` 与 `MECP.x` 的对应关系

驱动器在每一步对求解器生成如下形式的命令：

```text
./MECP.x astepN --algo <algorithm> --tde <tde> --tgmax <tgmax> --tgrms <tgrms> --tdxmax <tdxmax> --tdxrms <tdxrms> --stpmx <stpmx> --opt <opt_method> --hupd <hupd_method> --ndiis <ndiis> --pf_alpha <pf_alpha> --pf_sigma <pf_sigma> --pf_tstep <pf_tstep> --pf_tgrad <pf_tgrad>
```

因此，输入文件中的 `%Control` 是 `MECP.x` 参数的上层描述接口。

# 工作目录、状态文件与生成物

## 每一步的常见生成物

下表给出逐步工作目录中最常见的文件：

| 文件 | 类型 | 说明 |
| --- | --- | --- |
| `astepN.xyz` | 文本 | 第 `N` 步几何。 |
| `astepN_state1.<suffix>` | 文本 | state1 输入文件或脚本。 |
| `astepN_state2.<suffix>` | 文本 | state2 输入文件或脚本。 |
| `astepN_state1.<OUTPUT_SUFFIX>` | 文本 | state1 输出。 |
| `astepN_state2.<OUTPUT_SUFFIX>` | 文本 | state2 输出。 |
| `astepN.grad1` | 文本 | state1 能量与梯度。 |
| `astepN.grad2` | 文本 | state2 能量与梯度。 |
| `MECP.x` | 二进制 | 运行时编译得到的求解器。 |
| `MECP_temp.f` | 文本 | 按原子数改写后的 Fortran 源文件。 |
| `MECP.state` | 文本 | 重启状态文件。 |
| `MECP.diis` | 文本 | DIIS 历史文件，仅在 GDIIS / GEDIIS 模式下使用。 |
| `convg.tmp` | 文本 | 当前步收敛信息。 |
| `new.xyz` | 文本 | 未收敛时生成的下一步几何。 |
| `final.xyz` | 文本 | 收敛时生成的最终几何。 |
| `opt.trj` | 文本 | 轨迹文件。 |

## `.grad` 文件格式

`astepN.grad1` / `astepN.grad2` 当前格式如下：

1. 第 1 行：原子数；
2. 第 2 行：该态能量（Hartree）；
3. 后续每行：`Element  gx  gy  gz`。

当前结果链路要求这些梯度分量在写入 `.grad` 时对应 Hartree/Bohr 体系；求解器内部会统一换算为 Hartree/Angstrom。

## `convg.tmp`

`convg.tmp` 是驱动器和求解器之间的主要收敛状态文件。其内容由求解器在每一步重写，常见结构如下：

```text
NOT_CONVERGED
Step: 7
Energy_State_1:   ...
Energy_State_2:   ...
Energy_Gap:       ...
Convergence_Criteria:
Max_Gradient: ...
RMS_Gradient: ...
Max_Displacement: ...
RMS_Displacement: ...
Energy_Gap_Conv: ...
Parallel_Gradient_RMS: ...
Perpendicular_Gradient_RMS: ...
```

在 PF / Lagrange 模式下，还会追加：

```text
PF_Convergence_Criteria:
PF_Objective: ...
PF_Function1: ...
PF_Function2: ...
PF_Function3: ...
```

### 驱动器如何使用 `convg.tmp`

- 读取第一行判断是否为 `CONVERGED`；
- 从 `Step: N` 中恢复重启步号；
- 解析能量差、梯度与位移指标，生成 `opt.trj` 的注释行摘要。

## `opt.trj`

`opt.trj` 为 XYZ 轨迹文件。当前版本在追加每一帧时会改写第二行为收敛摘要，例如：

```text
E=0.000123 MaxF=0.001234 RMSF=0.000567 MaxD=0.003210 RMSD=0.001230
```

对于初始帧，在尚无 `convg.tmp` 时会写入：

```text
Convergence info not available
```

## 最终输出

收敛后，用户直接关心的输出通常为：

| 文件 | 说明 |
| --- | --- |
| `<geom_stem>_final.xyz` | 最终交叉点结构。 |
| `opt.trj` | 全程轨迹。 |

其中 `<geom_stem>_final.xyz` 的生成规则如下：

- 使用 `tmpdir` 运行时，程序将 `final.xyz` 复制回原目录并重命名；
- 直接在当前目录运行时，程序将 `final.xyz` 重命名为 `<geom_stem>_final.xyz`。

# 重启、临时目录与清理行为

## `tmpdir`

### `tmpdir=.`

- 在当前目录直接运行；
- 所有 `astep*`、输出日志、`.grad`、状态文件都写入当前目录。

### `tmpdir!=.`

- 进入指定目录运行；
- 若该目录已存在，程序会先将其整个删除，再创建新目录；
- 自动复制以下文件到工作目录：
  - `%Geom` 指定的几何文件；
  - `baneMECP.f`；
  - 选中的 `<prog>.conf`（若存在）。

除上述文件外，模板依赖的脚本或附属文件需要由模板脚本自行处理。

## `keeptmp`

### `keeptmp=true`

- 保留工作目录中的中间文件；
- 在 `tmpdir!=.` 时，计算结束后不删除临时目录。

### `keeptmp=false`

行为取决于当前运行位置：

#### 在临时目录中运行

程序会先将以下结果复制回原目录：

- `opt.trj`
- `<geom_stem>_final.xyz`（仅在收敛时）

随后删除整个 `tmpdir`。

#### 在当前目录中运行

程序执行清理操作。当前实现会删除的常见文件包括：

- `MECP.x`
- `MECP_temp.f`
- `*.log`
- `*.gms`
- `*.gjf`
- `*.chk`
- `new.xyz`
- `final.xyz`

在“非重启保留”清理路径中，还会删除：

- `astep*`
- `MECP.state`
- `convg.tmp`
- `opt.trj`

因此，在当前目录运行且 `keeptmp=false` 时，中间轨迹和状态文件不会保留。

## `restart`

### 判定条件

当同时满足以下条件时，程序按重启模式启动：

- `%Control restart=true`
- 工作目录中存在 `MECP.state` 或 `convg.tmp`

### 重启步号恢复

程序从 `convg.tmp` 中查找：

```text
Step: N
```

并从 `N+1` 步继续执行。

### 重启启动时的保留与清理

重启初始化时，程序保留以下对象：

- `MECP.state`
- `convg.tmp`
- `opt.trj`
- `astep*`

同时会清理常见输入/输出工件：

- `*.log`
- `*.gms`
- `*.gjf`
- `*.chk`
- `new.xyz`
- `final.xyz`

因此，重启依赖的是状态文件和步结构，而不是旧的日志或检查点文件。

## 新任务启动时的清理

当不进入重启模式时，程序会按“全新任务”清理路径删除旧文件，包括：

- `astep*`
- `MECP.state`
- `convg.tmp`
- `opt.trj`
- `MECP.x`
- `MECP_temp.f`
- 常见日志和输入后缀文件

# 示例

## 示例一：Gaussian 两态内置模式

以下示例使用 Gaussian 作为两态后端，在每一步生成两份 `.gjf` 并自动提取输出中的力：

```text
%Prog g16
%Geom ch2.xyz

%Control
maxcyc=50
tmpdir=tmp
keeptmp=true
restart=false
algorithm=harvey
opt_method=bfgs
hupd_method=bfgs
tde=5e-5
tgmax=7e-4
tgrms=5e-4
tdxmax=4e-3
tdxrms=2.5e-3
end

%InpTmplt1
%chk=state1.chk
#p B3LYP/def2SVP force [guess(read)]

State 1

0 1
[geom]


end

%InpTmplt2
%chk=state2.chk
#p B3LYP/def2SVP force [guess(read)]

State 2

0 3
[geom]


end

%GrpTmplt
state1=scf
state2=scf
end
```

配套 `g16.conf` 可写为：

```ini
[main]
env='source ~/.bashrc
module load g16'
INPUT_SUFFIX=gjf
OUTPUT_SUFFIX=log
RUN_CMD='g16 ${ACTUAL_INP} > ${ACTUAL_OUT}'

[scf]
E=SCF Done:\s+E\([^)]*\)\s*=\s*([-\d\.]+)
GRAD.Loacte=Forces (Hartrees/Bohr)
GRAD.NLineSkip=2
GRAD.TargetColumns=3,4,5
GRAD.EndBy=----
GRAD.Type=force
```

运行：

```bash
banemecp run.inp
```

## 示例二：external 模式 + 内置提取

以下示例使用脚本调用 GAMESS，脚本负责生成 `.gms` 输出，而驱动器负责从输出中提取两态梯度：

```text
%Prog gamess
%Geom ro.xyz

%Control
maxcyc=50
tmpdir=temp
keeptmp=true
restart=false
algorithm=bpupd
stpmx=0.05
end

%InpTmplt1
echo ">>>>> Start gamess calculation1"
cp ../gamess.sh .
bash gamess.sh [xyzfile] --iroot 1 --basename astep[nstep]_state1 >> mrsfout.tmp 2>&1
end

%InpTmplt2
echo ">>>>> Start gamess calculation2"
cp ../gamess.sh .
bash gamess.sh [xyzfile] --iroot 2 --basename astep[nstep]_state2 >> mrsfout.tmp 2>&1
end

%GrpTmplt
state1=mrsf1
state2=mrsf2
end
```

与之对应的 `gamess.conf` 可只保留 `OUTPUT_SUFFIX` 和状态段，而不写 `RUN_CMD`：

```ini
[main]
OUTPUT_SUFFIX=gms

[mrsf1]
E=1\s{2}A\s+(-?\d+\.?\d+)
GRAD.Loacte=GRADIENT OF THE ENERGY
GRAD.LoacteCount=2
GRAD.NLineSkip=3
GRAD.TargetColumns=3,4,5
GRAD.EndBy=
GRAD.Type=grad

[mrsf2]
E=2\s{2}A\s+(-?\d+\.?\d+)
GRAD.Loacte=GRADIENT OF THE ENERGY
GRAD.LoacteCount=2
GRAD.NLineSkip=3
GRAD.TargetColumns=3,4,5
GRAD.EndBy=
GRAD.Type=grad
```

## 示例三：external 模式 + 脚本直接写 `.grad`

以下示例中，脚本不再依赖 `%GrpTmplt`，而是在单个脚本中完成两态计算并直接写出 `.grad1` / `.grad2`：

```text
%Prog external
%Geom mol.xyz

%Control
maxcyc=30
tmpdir=tmp
keeptmp=true
restart=false
end

%InpTmplt1
python make_two_state_grad.py [xyzfile] [name]
end
```

在该模式下，`make_two_state_grad.py` 应保证在每一步结束后生成：

- `astepN.grad1`
- `astepN.grad2`

并且无需提供 `%InpTmplt2` 与 `%GrpTmplt`。

# 当前行为与限制

## 输入文件与路径

- 相对路径以程序启动时的当前工作目录为基准，不以输入文件所在目录为基准。
- `%Geom` 当前只支持标准 XYZ。
- 输入文件本身不定义稳定的注释语法；控制段建议不要在 `key=value` 行尾追加注释。
- 段结束符必须写为独立的 `end`。

## 程序模式与两态要求

- baneMECP 是两态 MECP 驱动程序，不支持三态或更多态的联合优化。
- 在内置模式中，程序总是需要 state1 与 state2 两套结果。
- 当使用 `%GrpTmplt` 启用内置提取时，程序总是提取 state1 与 state2；该工作流要求两态输出都存在。
- `%Prog external` 仅在模板脚本自行生成 `.grad1` / `.grad2` 时可独立使用；若要在 external 模式下继续使用内置提取，应提供可解析的 `<prog>.conf`。

## 配置文件

- `[main]` 中 `env`、`INPUT_SUFFIX`、`OUTPUT_SUFFIX`、`RUN_CMD` 使用固定键名。
- 当前驱动器实际读取的提取键名为 `E.LoacteCount`、`GRAD.Loacte`、`GRAD.LoacteCount`；其他拼写不会进入当前结果提取流程。
- `GRAD.unit` 当前不驱动单位切换逻辑；源输出应按 Hartree/Bohr 组织。
- `%GrpTmplt` 只读取前两行，且按 `state1` 在前、`state2` 在后的顺序解析。
- `%GrpTmplt` 所引用的 section 名需要与配置文件完全一致，并且建议仅使用字母、数字和下划线。

## 提取与梯度块读取

- 能量提取只使用正则表达式的第一个捕获组。
- 梯度块提取在定位并跳过指定行数后，按原子数固定读取梯度行数。
- `GRAD.EndBy` 用于检查是否提前遇到结束标记或空行，不作为可变长度扫描终止条件。
- 若梯度块中夹杂额外说明行、分页行或重复表头，且这些行落在预期原子数据范围内，当前提取流程会失败。

## shell 与运行环境

- external 模式固定以 `bash` 执行脚本。
- 内置模式中的 `env` 与 `RUN_CMD` 通过系统默认 shell 执行，命令内容应与该 shell 兼容，或显式调用所需 shell。
- 程序在运行时需要 `gfortran`；缺少该编译器时无法生成 `MECP.x`。

## 临时目录与清理

- 当 `tmpdir` 不为 `.` 时，若该目录已存在，启动时会被整个删除并重建。
- 当在当前目录运行且 `keeptmp=false` 时，中间文件、状态文件与 `opt.trj` 会被清理。
- 当在临时目录运行且 `keeptmp=false` 时，只复制回 `opt.trj` 和最终结构；其他中间文件不保留。

# 附录：完整控制段示例

```text
%Control
maxcyc=100
tmpdir=tmp
keeptmp=true
debug=false
restart=false

algorithm=harvey
opt_method=bfgs
ndiis=4
hupd_method=bfgs

tde=5e-5
tgmax=7e-4
tgrms=5e-4
tdxmax=4e-3
tdxrms=2.5e-3
stpmx=0.1

pf_alpha=0.02
pf_sigma=3.5
pf_tstep=1e-6
pf_tgrad=5e-3
end
```

# 附录：典型输出文件示意

```text
case01/
  run.inp
  mol.xyz
  mol_final.xyz
  opt.trj
  tmp/
    astep1.xyz
    astep1_state1.gjf
    astep1_state1.log
    astep1_state2.gjf
    astep1_state2.log
    astep1.grad1
    astep1.grad2
    astep2.xyz
    ...
    MECP.x
    MECP_temp.f
    MECP.state
    MECP.diis
    convg.tmp
```
