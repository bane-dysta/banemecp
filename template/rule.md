- `%Prog`:
定义输入程序。没有配置文件或有配置文件却无执行命令(RUN_CMD)的全部视作external类。用于加载`%InpTmplt`执行信息。


- `%Geom`:
输入初始结构，要求标准xyz文件。

- `%Control`
  - `maxcyc(参数：数值)`: 最大迭代圈数，默认100 
  - `tmpdir(参数：路径)`: 临时文件目录
  - `keeptmp(参数：true/false)`: 是否保留临时文件
  - `restart(参数：true/false)`: 是否重启


- `%InpTmplt1`、`%InpTmplt2`
  第1、2个态的输入模板。内置程序的启动命令需要在配置文件给出，external类则直接eval执行从`%InpTmplt`开始到end内的命令。
  为实现迭代，本模块使用各种占位符:
  - `[name]` - 替换为当前将要传递给MECP求解器的输入文件名 (step1,2,...,n,无后缀);
  - `[xyzfile]` - 替换为当前结构的xyz文件名 (step1,2,...,n.xyz,等同于[name].xyz);
  - `[geom]` - 替换为当前结构坐标;
  - `[nstep]` - 替换为当前步数
  - `[tmpdir]` - 替换为tmpdir指定的临时文件存储路径
  - `[其他占位符]`：视作step1忽略，step2开始替换为方框内的内容，用于从第二步开始启用guess=read等

  注意：若external模式下使用tmpdir，你可能需要将external中涉及到的脚本复制到tmpdir。external模式下shell的启动位置是tmpdir内。

- `%GrpTmplt`
  使用内置梯度/能量提取逻辑，分两行，第一行为state1类型，第二行为state2类型，要求值在对应配置文件中有定义;
  external类程序如果定义了提取逻辑，可以使用`%GrpTmplt`，否则需要`%InpTmplt`中执行的命令自己产生.grad1、.grad2文件。其格式等同于xyz文件，但坐标替换为梯度值(Hartrees/Bohr)，第二行写入能量(Hartree)。

