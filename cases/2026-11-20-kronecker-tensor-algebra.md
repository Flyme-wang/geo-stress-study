# 案例拆解：《连续介质力学（二）：Kronecker 符号与二阶张量运算》

> 来源：微信公众号「学以远志」，2026-11-20 20:15
>
> 篇幅：一篇指标记法长文，11 节、4 道练习；公式以图片嵌入，共 173 张，无独立技术插图
>
> 对应知识图谱节点：`C7`
>
> 用途：应力张量之前的数学预备；把"指标、基底、映射、双点积"一次讲清
>
> 本库正文：[M0-04 Kronecker 符号与二阶张量](../docs/00-notation/kronecker-and-tensor-algebra.md)

---

## 0. 结论先行

这篇推文的定位很明确：**在讲应力之前，先把张量代数的"语法"补齐**。它从三维正交单位基底出发，依次讲 Kronecker 符号、矢量张量积、二阶张量的九个基张量、单位张量、乘法、转置与对称/反对称分解、迹与双点积，最后用一道 $\mathbf{S}:\mathbf{W}=0$ 的例题和四道练习收尾。链条完整，顺序合理，适合作为应力张量笔记的前置阅读。

逐项核对后，原文的数学结论没有发现错误。 $\delta\_{ij}a\_j=a\_i$、 $(\mathbf{a}\otimes\mathbf{b})\mathbf{c}=\mathbf{a}(\mathbf{b}\cdot\mathbf{c})$、 $C\_{ij}=A\_{ik}B\_{kj}$、 $\mathbf{A}=\mathbf{S}+\mathbf{W}$、 $\mathbf{A}:\mathbf{B}=A\_{ij}B\_{ij}$ 都与独立推导一致。它把"张量积不交换""双点积有两种配对"这些容易含糊的地方单独点出来，这一点比多数教材处理得好。

本库补强的是四件事：**指标规则可执行**（重复/自由指标逐条判定）、**矩阵与张量分清**（矩阵只是基底下的分量表）、**双点积配对写死**（同位 vs 交换，配数值反例）、**全部数字可复算**（参考实现与校验脚本）。另外，原文的公式图片丢失了若干单字母符号，本库按上下文补全，记录在 §3.1。

适用范围：三维实矢量、实二阶张量、正交单位基底、列矢量约定。原文与本文都不涉及高阶张量、曲线坐标、协变导数与张量场求导；这些内容需要另开页面。

---

## 1. 原文结构与论证链

| 阶段 | 原文讲什么 | 本库位置 |
|---|---|---|
| 一、适用范围 | 三维直角坐标、正交单位基底、求和约定 | [M0-04 §1.1](../docs/00-notation/kronecker-and-tensor-algebra.md) |
| 二、Kronecker 符号 | 定义、指标替换、三种缩并 | [M0-04 §2.1](../docs/00-notation/kronecker-and-tensor-algebra.md) |
| 三、矢量的张量积 | 外积矩阵、作用规则、非交换性 | [M0-04 §2.2](../docs/00-notation/kronecker-and-tensor-algebra.md) |
| 四、基底表示 | 九个基张量、基张量的作用、分量运算 | [M0-04 §2.3](../docs/00-notation/kronecker-and-tensor-algebra.md) |
| 五、单位张量 | 分量、作用、与任意张量的乘法 | [M0-04 §2.4](../docs/00-notation/kronecker-and-tensor-algebra.md) |
| 六、二阶张量乘法 | 映射复合、行乘列、非交换性 | [M0-04 §2.5](../docs/00-notation/kronecker-and-tensor-algebra.md) |
| 七、转置与分解 | 转置规则、对称/反对称、唯一分解 | [M0-04 §2.6](../docs/00-notation/kronecker-and-tensor-algebra.md) |
| 八、迹与双点积 | 迹、基张量内积、Frobenius 范数、两种缩并 | [M0-04 §2.7](../docs/00-notation/kronecker-and-tensor-algebra.md) |
| 九、例题 | $\mathbf{S}:\mathbf{W}=0$ | [M0-04 §2.8](../docs/00-notation/kronecker-and-tensor-algebra.md) |
| 十、练习 | 化简、张量积、两种缩并、分解 | [M0-04 §5](../docs/00-notation/kronecker-and-tensor-algebra.md) |

这张表也回答了"读原文还是读本库"：原文适合第一遍建立框架，公式排列紧凑；要查某个指标为什么被筛掉、要写代码复算、要准备课堂练习，用本库对应小节与校验脚本。

---

## 2. 原文讲对的四处

### 2.1 先定义基底，再谈分量

原文没有直接写 $A\_{ij}$，而是先固定三维正交单位基底，说明 $\mathbf{e}\_i\cdot\mathbf{e}\_j=\delta\_{ij}$，再引入分量。这个顺序是对的：分量的意义依赖基底；不先固定基底，"第 $i$ 行第 $j$ 列"就没有几何解释。本库进一步把分量定义为 $A\_{ij}=\mathbf{e}\_i\cdot(\mathbf{A}\mathbf{e}\_j)$，让"输入沿哪个基矢量、输出取哪个分量"一目了然。

### 2.2 指标替换被讲成了动作，而不是公式

原文用" $\delta$ 在求和中筛选出对应分量"来解释 $\delta\_{ij}a\_j=a\_i$，并给出三个 $i$ 的展开。这比只写一句" $\delta$ 是单位矩阵"更有用：学生能看到"固定 $i$、对 $j$ 求和、只有 $i=j$ 的项留下"这三步。本库把这三步画成图 51(b)，并把 $\delta\_{11}=1$ 与 $\delta\_{ii}=3$ 的差别单独列出。

### 2.3 张量积的作用规则没有跳步

原文先写 $(\mathbf{a}\otimes\mathbf{b})\_{ij}=a\_i b\_j$，再写 $(\mathbf{a}\otimes\mathbf{b})\mathbf{c}=\mathbf{a}(\mathbf{b}\cdot\mathbf{c})$，并用分量推导把"提出公共因子 $a\_i$"的过程展开。这个推导是理解"张量积 = 先投影、再定向输出"的关键。本库补了一张几何图（图 52(b)），把标量 $\mathbf{b}\cdot\mathbf{c}$ 与输出方向 $\mathbf{a}$ 分开画。

### 2.4 双点积的两种配对被明确区分

原文没有把 $A\_{ij}B\_{ij}$ 与 $A\_{ij}B\_{ji}$ 混为一谈，而是并列写出，并用一个只有单个非零元的例子说明两者不同。这是本文最有价值的一处：很多教材只写 $\mathbf{A}:\mathbf{B}=A\_{ij}B\_{ij}$，读者遇到 $\mathrm{tr}(\mathbf{A}\mathbf{B})$ 时就不知道对应哪种配对。本库把这个反例扩展成图 55，并补上 $\mathrm{tr}(\mathbf{A}^{\mathsf T}\mathbf{B})$ 与 $\mathrm{tr}(\mathbf{A}\mathbf{B})$ 的对照。

---

## 3. 本库补强与校核的六处

### 3.1 公式图片丢失单字母符号：源格式问题

原文的公式以图片嵌入，HTML 的 `alt` 文本保存了 LaTeX，但正文中若干单字母符号在图片之外的位置丢失。例如原文出现"其中，□ 是哑指标""第 □ 行、第 □ 列的元素"这样的空位（□ 表示原文缺失处）；从上下文看，空位分别是 $i$ 、 $j$ 等。本库在转写时按上下文补全，并在所有公式中显式写出指标。

这不是原文的数学错误，而是公众号排版把行内符号拆成图片后留下的格式问题。它提醒我们：**公式图片不适合作为唯一载体**；一旦 `alt` 或正文丢字，读者无法从上下文恢复。本库因此坚持用 LaTeX 源文，并对每个公式做数值核对。

### 3.2 重复指标与自由指标：一条可执行的判定规则

原文讲了求和约定，但没有把"同一项里同一字母最多出现两次"写成硬规则。本库补上这条规则，并给出三种情形的对照： $a\_i b\_i$ 是标量（ $i$ 重复）， $a\_i$ 是矢量分量（ $i$ 自由）， $A\_{ij}$ 是二阶张量分量（ $i,j$ 都自由）。出现三次指标时，表达式有歧义，必须改写或加括号。

这条规则在应力笔记里会反复用到： $\sigma\_{ij}n\_j$ 是一个自由指标 $i$ 的矢量； $\sigma\_{ij}n\_i n\_j$ 没有自由指标，是标量； $\sigma\_{ii}$ 是迹。看到表达式先数指标，比背公式可靠。

### 3.3 矩阵不是张量：换基时什么变、什么不变

原文为了计算方便，基本把二阶张量等同于 $3\times3$ 矩阵。对固定基底的运算是安全的，但学生容易由此形成"张量就是矩阵"的印象。本库明确区分：张量是线性映射，矩阵是某组基下的分量表；换基时矩阵按相似变换改变，映射本身不变。

由此可以解释两个常见困惑：为什么转置在某些基下看起来像"交换行列"，但它其实是交换输入/输出的配对；为什么迹与 Frobenius 范数在正交变换下不变，而单个分量会变。前者是映射的性质，后者是分量的性质。

### 3.4 双点积的配对必须写死

原文并列了两种缩并，但没有给出一般的迹表达式。本库补上

$$
\mathbf{A}:\mathbf{B}=A_{ij}B_{ij}=\mathrm{tr}(\mathbf{A}^{\mathsf T}\mathbf{B}),
\qquad
A_{ij}B_{ji}=\mathrm{tr}(\mathbf{A}\mathbf{B}).
$$

并用 $\mathbf{A}=\mathbf{B}=\mathbf{e}\_1\otimes\mathbf{e}\_2$ 验证：第一种得 1，第二种得 0。这个反例很小，但足以说明"都是标量"不等于"同一个标量"。在应力功共轭、能量范数与切线刚度里，配错一次就会让能率差一个转置。

### 3.5 正交变换下的不变量：补一条可证伪的检验

原文说 Frobenius 范数在正交基变换下不变，但没有给数值。本库取转轴 $[1,2,-1]$、转角 37 度，构造旋转矩阵 $\mathbf{Q}$，验证：

| 检验 | 结果 |
|---|---|
| $\det\mathbf{Q}$ | 1 |
| $\mathbf{Q}^{\mathsf T}\mathbf{Q}-\mathbf{I}$ 的最大元素 | $7.99\times10^{-17}$ |
| $\mathrm{tr}(\mathbf{Q}\mathbf{A}\mathbf{Q}^{\mathsf T})-\mathrm{tr}\mathbf{A}$ | 0 |
| $\Vert\mathbf{Q}\mathbf{A}\mathbf{Q}^{\mathsf T}\Vert\_{\mathrm F}-\Vert\mathbf{A}\Vert\_{\mathrm F}$ | 0 |

这组数字的作用是"证伪"：如果有人把范数写成"分量平方和"却忘了开方，或者把迹当成 $\mathbf{A}:\mathbf{A}$，这里的差值会立刻暴露。

### 3.6 课堂练习：从化简到反例

原文的四道练习覆盖了化简、张量积、两种缩并与分解。本库保留这四道，另加一道判断题：若 $\mathbf{A}:\mathbf{B}=0$，是否必有 $\mathbf{A}=\mathbf{0}$ 或 $\mathbf{B}=\mathbf{0}$？答案是否定的： $\mathbf{e}\_1\otimes\mathbf{e}\_1$ 与 $\mathbf{e}\_2\otimes\mathbf{e}\_2$ 双点积为零，但两者都非零。这道题把"双点积是内积"与"内积为零意味着正交，不意味着零元"区分开，为后面讲应力空间的正交分解做铺垫。

---

## 4. 数值对照

全部数字由 [`tools/verify_tensor_notes.py`](../tools/verify_tensor_notes.py) 复算，末尾打印 `ALL CHECKS PASSED`。下表只列与本案例直接相关的几组。

| 检验 | 输入 | 结果 |
|---|---|---|
| 指标替换 | $\mathbf{a}=[1,2,3]$ | $\delta\_{ij}a\_j=[1,2,3]$ |
| 两种缩并 | $\delta\_{ii}$、 $\delta\_{ij}\delta\_{ij}$ | 3、3 |
| 张量积作用 | $\mathbf{a}=[1,2,3]$、 $\mathbf{b}=[4,5,6]$、 $\mathbf{c}=[-1,2,1]$ | $\mathbf{b}\cdot\mathbf{c}=12$，输出 $[12,24,36]$ |
| 非交换性 | 同上 | $\max\vert\mathbf{a}\otimes\mathbf{b}-\mathbf{b}\otimes\mathbf{a}\vert=6$ |
| 基张量作用 | $\mathbf{v}=[2,-3,5]$ | $(\mathbf{e}\_1\otimes\mathbf{e}\_2)\mathbf{v}=[-3,0,0]$ |
| 乘法非交换 | 文中的 $\mathbf{A}$、 $\mathbf{B}$ | $\max\vert\mathbf{A}\mathbf{B}-\mathbf{B}\mathbf{A}\vert=2$ |
| 对称/反对称 | 文中 $\mathbf{A}$ | $\mathbf{S}:\mathbf{W}=0$，随机矩阵分解误差 $2.22\times10^{-16}$ |
| 两种缩并反例 | $\mathbf{A}=\mathbf{B}=\mathbf{e}\_1\otimes\mathbf{e}\_2$ | $A\_{ij}B\_{ij}=1$， $A\_{ij}B\_{ji}=0$ |
| 迹与范数 | 文中 $3\times3$ 矩阵 | $\mathrm{tr}=8$， $A\_{ij}A\_{ij}=30$， $A\_{ij}A\_{ji}=26$ |
| 旋转不变性 | 轴 $[1,2,-1]$、37 度 | $\det\mathbf{Q}=1$，正交性误差 $7.99\times10^{-17}$ |

这些数字的用途不是"再算一遍"，而是给每个公式配一条独立检查：矩阵乘法、分量求和、几何投影、旋转不变性，四条路径至少两条对上，才认为公式写对了。

---

## 5. 图件与版权登记

### 5.1 原文图件

原文没有独立的技术插图。它的"图"全部是公式图片：`#js_content` 中共 421 个 DOM 元素，其中 173 个是公式图片，每个图片的 `alt` 保存了对应的 LaTeX；`block` 标记区分行内公式与独立公式。文章没有实验照片、示意图或数据图。

本库没有复制任何公式图片。正文中的公式全部按 `alt` 的 LaTeX 重新排版，并做了两件事：补全丢失的单字母符号；对每个可计算的结论跑数值校验。

### 5.2 本库新增图件

| 图 | 文件 | 内容 |
|---|---|---|
| 图 51 | `assets/figures/fig51-kronecker-index.png` | $\delta\_{ij}$ 的分量矩阵、指标替换、两种缩并 |
| 图 52 | `assets/figures/fig52-tensor-product.png` | 外积矩阵、投影作用、非交换差矩阵 |
| 图 53 | `assets/figures/fig53-basis-tensors.png` | 九个基张量、基张量作用、分量运算 |
| 图 54 | `assets/figures/fig54-symmetric-decomposition.png` | $\mathbf{A}=\mathbf{S}+\mathbf{W}$ 的热图分解 |
| 图 55 | `assets/figures/fig55-double-dot.png` | 两种缩并的配对与 Frobenius 范数 |

五张图均由 [`tools/make_figures.py`](../tools/make_figures.py) 生成，未使用原文图片。图 52(b) 的几何构图为教学示意，矢量的长度按显示比例缩短，不表示原文的数值比例。

### 5.3 版权边界

原文的文字、公式图片与练习表述版权归原作者。本库只吸收其知识框架与题目类型，正文独立撰写；公式按数学内容重新排版，数值与图件全部重做。原文链接仅登记在案例笔记与正文来源中，不嵌入原文图片。若需在课堂使用原文，请通过公众号原文链接访问。

---

## 6. 课堂用法与图谱映射

### 6.1 两学时安排

| 时段 | 内容 | 对应图 |
|---|---|---|
| 0–15 min | 指标、自由/哑指标、 $\delta$ 的筛选与计数 | 图 51 |
| 15–35 min | 张量积、外积矩阵、先投影再输出 | 图 52 |
| 35–55 min | 九个基张量、基底表示、分量运算 | 图 53 |
| 55–75 min | 乘法、转置、对称/反对称分解 | 图 54 |
| 75–95 min | 迹、双点积、两种配对、Frobenius 范数 | 图 55 |
| 95–120 min | 四道练习 + 判断题，现场用脚本复算 | — |

课堂建议：先让学生用手算 $\delta\_{ij}a\_j$ 与 $\mathbf{a}\otimes\mathbf{b}$，再打开脚本对照。图 52(b) 与图 53(b) 适合投屏讲"读出分量、定向输出"的动作；图 54 与图 55 适合讲"分解"与"配对"的区别。

### 6.2 图谱映射

- 模块：`M0` 数学与约定
- 主题节点：`M0-04` Kronecker 符号与二阶张量
- 案例节点：`C7` 案例：Kronecker 符号与二阶张量
- 前置边：`M0-01`（张量与指标记法）→ `M0-04`
- 包含边：`M0` → `M0-04`
- 案例边：`C7` → `M0-04`

新增节点后，图谱从 64 个节点、106 条边变为 66 个节点、109 条边。页面数据由 [`tools/build_graph_data.py`](../tools/build_graph_data.py) 重新生成。

---

## 7. 来源

- 原文：公众号「学以远志」，《连续介质力学（二）：Kronecker 符号与二阶张量运算》，2026-11-20 20:15。链接：<https://mp.weixin.qq.com/s/aqgMf6rpw3EGnYaSFg5P3Q>。
- 原文结构：11 节（适用范围；Kronecker 符号；矢量张量积；基底表示；单位张量；乘法；转置、对称与反对称；迹与双点积；例题；练习；练习答案），4 道练习。
- 本库正文：[`docs/00-notation/kronecker-and-tensor-algebra.md`](../docs/00-notation/kronecker-and-tensor-algebra.md)
- 参考实现：[`tools/tensor_lib.py`](../tools/tensor_lib.py)
- 校验脚本：[`tools/verify_tensor_notes.py`](../tools/verify_tensor_notes.py)
- 图件脚本：[`tools/make_figures.py`](../tools/make_figures.py)
