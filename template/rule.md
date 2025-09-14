# baneMECP 输入格式规则

## 基本指令

- `%Prog`:
定义输入程序。没有配置文件或有配置文件却无执行命令(RUN_CMD)的全部视作external类。用于加载`%InpTmplt`执行信息。

- `%Geom`:
输入初始结构，要求标准xyz文件。

## 控制参数

- `%Control`
  
  **基本控制参数：**
  - `maxcyc(参数：数值)`: 最大迭代圈数，默认50
  - `tmpdir(参数：路径)`: 临时文件目录，默认当前目录(".")
  - `keeptmp(参数：true/false)`: 是否保留临时文件，默认false
  - `debug(参数：true/false)`: 是否输出调试信息，默认false
  - `restart(参数：true/false)`: 是否重启，默认false
  
  **收敛判据参数：**
  - `tde(参数：数值)`: 能量间隙收敛阈值，默认5e-5 (Hartree)
  - `tgmax(参数：数值)`: 最大梯度收敛阈值，默认7e-4 (Hartree/Angstrom)
  - `tgrms(参数：数值)`: RMS梯度收敛阈值，默认5e-4 (Hartree/Angstrom)
  - `tdxmax(参数：数值)`: 最大位移收敛阈值，默认4e-3 (Angstrom)
  - `tdxrms(参数：数值)`: RMS位移收敛阈值，默认2.5e-3 (Angstrom)
  
  **优化控制参数：**
  - `stpmx(参数：数值)`: 置信半径/最大步长，默认0.1 (Angstrom)
    - 对于振荡系统，建议减小至0.02-0.05以提高稳定性
    - 该参数控制每步的最大几何变化，影响优化的稳定性和速度

## 输入模板

- `%InpTmplt1`、`%InpTmplt2`
  第1、2个态的输入模板。内置程序的启动命令需要在配置文件给出，external类则直接eval执行从`%InpTmplt`开始到end内的命令。
  
  **占位符系统：**
  - `[name]` - 替换为当前将要传递给MECP求解器的输入文件名 (step1,2,...,n,无后缀)
  - `[xyzfile]` - 替换为当前结构的xyz文件名 (step1,2,...,n.xyz,等同于[name].xyz)
  - `[geom]` - 替换为当前结构坐标
  - `[nstep]` - 替换为当前步数
  - `[tmpdir]` - 替换为tmpdir指定的临时文件存储路径
  - `[其他占位符]`：视作step1忽略，step2开始替换为方框内的内容，用于从第二步开始启用guess=read等

  **注意事项：**
  - `%InpTmplt2`为可选项，如果未提供则只进行单态计算
  - external模式下使用tmpdir时，需要将相关脚本复制到tmpdir
  - external模式下shell的启动位置是tmpdir内

## 梯度提取

- `%GrpTmplt`
  使用内置梯度/能量提取逻辑，分两行：
  - 第一行：state1类型
  - 第二行：state2类型
  
  要求值在对应配置文件中有定义。external类程序如果定义了提取逻辑，可以使用`%GrpTmplt`，否则需要在`%InpTmplt`中执行的命令自己产生.grad1、.grad2文件。
  
  **梯度文件格式：**
  - 格式等同于xyz文件，但坐标替换为梯度值(Hartrees/Bohr)
  - 第二行写入能量(Hartree)
  - 程序会自动处理单位转换(Hartree/Bohr → Hartree/Angstrom)

## 示例配置

```
%Prog gamess
%Geom molecule.xyz
%Control
maxcyc=100
tmpdir=./tmp
keeptmp=false
restart=false
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

