# baneMECP 输入文件规则

输入文件由若干指令段组成；段名以`%`开头。

## 基本指令

- `%Prog <prog>`

  选择量子化学程序/驱动方式，并用于定位配置文件`<prog>.conf`。

  - `prog=external`：强制external模式；配置文件可选（仅在使用`%GrpTmplt`内置提取时才需要）。
  - `prog!=external`：必须能找到`<prog>.conf`；若找不到则报错退出。
  - 若`<prog>.conf`存在但`[main]`段缺少`RUN_CMD`，将自动切换为external模式。

- `%Geom <geom.xyz>`

  初始几何结构（标准xyz文件）。

## %Control 控制参数

`%Control`段使用`key=value`逐行给出，以单独一行`end`结束。

布尔值仅接受`true/false`（小写）。

### 基本控制

- `maxcyc=<int>`：最大迭代步数。默认`50`。
- `tmpdir=<path>`：计算临时目录。默认`.`。
  - `tmpdir=.`：在当前目录运行。
  - `tmpdir!=.`：进入该目录运行（若目录已存在，会被删除并重新创建）。
- `keeptmp=true/false`：是否保留临时目录/中间文件。默认`false`。
- `debug=true/false`：输出调试信息。默认`false`。
- `restart=true/false`：是否尝试从已有状态继续迭代。默认`false`。
  - `restart=true`且工作目录存在`MECP.state`或`convg.tmp`时，将从`convg.tmp`中读取`Step: N`并从`N+1`继续；否则从头开始。

### 收敛判据（传递给MECP.x）

- `tde=<float>`：能量差收敛阈值（Hartree）。默认`5e-5`。
- `tgmax=<float>`：最大梯度收敛阈值（Hartree/Angstrom）。默认`7e-4`。
- `tgrms=<float>`：RMS梯度收敛阈值（Hartree/Angstrom）。默认`5e-4`。
- `tdxmax=<float>`：最大位移收敛阈值（Angstrom）。默认`4e-3`。
- `tdxrms=<float>`：RMS位移收敛阈值（Angstrom）。默认`2.5e-3`。
- `stpmx=<float>`：信赖半径/最大步长（Angstrom）。默认`0.1`。

### MECP 算法选择

- `algorithm=<name>`：MECP求解算法。默认`harvey`。
  - `harvey`
  - `bpupd`
  - `lagrange`（别名：`lc`，`pf`）

### 优化加速与Hessian更新

- `opt_method=<name>`：步进/加速方法。默认`bfgs`。
  - `bfgs`
  - `gdiis`
  - `gediis`
  - 兼容别名：`opt`、`optmethod`、`opt_method`、`diis_method`。
- `ndiis=<int>`：DIIS子空间维数（历史长度）。默认`4`。
  - 兼容别名：`diis_size`、`diissize`。
- `hupd_method=<name>`：Hessian更新方式。默认`bfgs`。
  - `bfgs`
  - `psb`
  - `sr1`
  - `bofill`
  - 兼容别名：`hupd`、`hupd_method`、`hess_update`、`hessupd`。

### Lagrange / 惩罚函数参数

仅在`algorithm=lagrange`或`algorithm=pf`时生效。

- `pf_alpha=<float>`：alpha（Hartree）。默认`0.02`。
- `pf_sigma=<float>`：sigma（1/Hartree）。默认`3.50`。
- `pf_tstep=<float>`：罚函数目标值变化阈值。默认`1e-6`。
- `pf_tgrad=<float>`：罚函数梯度阈值。默认`5e-3`。
- `pf_thresh=<tstep,tgrad>`：同时设置`pf_tstep`与`pf_tgrad`（支持`1e-6,0.005`或`[1e-6, 0.005]`形式）。
  - 兼容键：`pf_thresh_step`、`pf_thresh_grad`。

## %InpTmplt 输入模板

- `%InpTmplt1`：必填。第1态输入/脚本模板。
- `%InpTmplt2`：可选。第2态输入/脚本模板。
- 每段以单独一行`end`结束。

### 模板执行方式

- **内置模式**（`<prog>.conf`存在且`[main]`含`RUN_CMD`）

  渲染后生成：
  - `astepN_state1.<INPUT_SUFFIX>`
  - `astepN_state2.<INPUT_SUFFIX>`（需要`%InpTmplt2`）

  并按`RUN_CMD`分别运行两态计算。

- **external模式**（`%Prog external`或配置无`RUN_CMD`）

  渲染后生成并执行：
  - `astepN_state1.sh`
  - `astepN_state2.sh`（需要`%InpTmplt2`）

### 占位符

在渲染模板时替换：

- `[name]`：当前步基名，形如`astep1`、`astep2`、…（无后缀）。
- `[xyzfile]`：当前步坐标文件名，形如`astep1.xyz`、`astep2.xyz`、…。
- `[geom]`：当前步几何坐标（xyz正文：元素符号 + x y z）。
- `[nstep]`：当前步号（1,2,...）。
- `[tmpdir]`：`tmpdir`参数的原始字符串值。

条件占位符（用于“从第二步开始生效”的内容）：

- 任意形如`[anything]`的片段：
  - 在第1步会被整体删除
  - 从第2步开始替换为`anything`（去掉方括号）

### %InpTmplt2 省略规则

MECP求解器需要同时存在`astepN.grad1`与`astepN.grad2`。

因此通常应同时提供`%InpTmplt1`与`%InpTmplt2`。仅在external模式且由`%InpTmplt1`自行完成两态计算并直接生成`.grad1/.grad2`时，`%InpTmplt2`才可省略。

## %GrpTmplt 梯度/能量提取选择

`%GrpTmplt`用于启用内置提取逻辑，以单独一行`end`结束。格式固定为两行：

- `state1=<type>`
- `state2=<type>`

其中`<type>`为配置文件`<prog>.conf`中的section名称。

- 提供`%GrpTmplt`：程序从输出文件`astepN_stateX.<OUTPUT_SUFFIX>`中按配置规则提取能量与梯度，并生成`astepN.grad1/.grad2`。
- 不提供`%GrpTmplt`：仅在external模式下可用，此时要求模板脚本自行生成`astepN.grad1`与`astepN.grad2`。

### .grad 文件格式

- 第1行：原子数
- 第2行：能量（Hartree）
- 后续每行：`Element  gx  gy  gz`（Hartree/Bohr）

程序会将梯度单位自动换算为Hartree/Angstrom。

## 示例

```
%Prog gamess
%Geom molecule.xyz
%Control
maxcyc=100
tmpdir=./tmp
keeptmp=false
restart=false
algorithm=harvey
opt_method=bfgs
hupd_method=bfgs
ndiis=4
tde=5e-5
tgmax=7e-4
tgrms=5e-4
tdxmax=4e-3
tdxrms=2.5e-3
stpmx=0.1
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
cp opt.trj ..
end

%GrpTmplt
state1=mrsf1
state2=mrsf2
end
```

