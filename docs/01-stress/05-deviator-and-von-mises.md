# 偏应力与等效应力：为什么三个主应力能压成一个数

> 对应知识图谱节点：`M1-05`
> 前置：[主应力](03-principal-stress.md)、[莫尔圆](04-mohr-circle.md)
> 后续：[屈服准则与岩石强度](../03-constitutive/README.md)（提纲）

---

## 0. 结论先行

1. 任何应力张量都能唯一分解为**体积部分**（球张量）与**形状改变部分**（偏张量）：

$$
\sigma_{ij}=\underbrace{p\,\delta_{ij}}_{\text{球张量}}+\underbrace{s_{ij}}_{\text{偏张量}},
\qquad p=\frac{\sigma_{kk}}{3}=\frac{\sigma_1+\sigma_2+\sigma_3}{3}
$$

2. 球张量只改变体积，不改变形状；**引起形状改变（畸变）的是偏张量**。
3. 偏张量的第二不变量

$$
J_2=\tfrac12 s_{ij}s_{ij}
=\frac{1}{6}\left[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2\right]
$$

是一个与坐标系无关的标量，度量"畸变有多强"。

4. **von Mises 等效应力**：

$$
\sigma_{\mathrm v}=\sqrt{3J_2}
=\sqrt{\tfrac12\left[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2\right]}
$$

5. 它是**一个假想的单轴应力值**，其畸变能/畸变程度与真实三维状态相同。
   它可以在单轴拉伸试验的屈服强度上做比较，但它**不是任何实际截面上的正应力**。

---

## 1. 前提与符号

- 各向同性材料；关注"屈服"这类由畸变驱动的过程。
- $p$：平均正应力（**不是**孔隙压力）。地质文献中的孔隙压力记作 $p_p$ 或 $p_f$，务必区分。
- $s_{ij}$：偏应力张量，$s_{kk}=0$（迹为零）。
- $J_2$：偏应力第二不变量；工程上也记作 $q$、$\sigma_{\mathrm{eq}}$、$\sigma_{\mathrm{vM}}$。

---

## 2. 推导

### 2.1 第一步：为什么可以这样拆

先证明球张量与偏张量正交（在张量内积意义下）：

$$
p\,\delta_{ij}\,s_{ij}=p\,s_{ii}=p\cdot 0=0
$$

**依据**：$\delta_{ij}s_{ij}=s_{ii}$，而 $s_{ii}=\sigma_{ii}-3p=3p-3p=0$。

正交意味着：**两部分互不影响**，可以分别讨论"体积效应"和"形状效应"。
这是分解有意义的根本原因（否则拆了也不能独立分析）。

### 2.2 第二步：$J_2$ 的三种写法

**（a）定义式**

$$
J_2=\tfrac12 s_{ij}s_{ij}
$$

**（b）用应力本身表示。** 把 $s_{ij}=\sigma_{ij}-p\delta_{ij}$ 代入：

$$
s_{ij}s_{ij}
=(\sigma_{ij}-p\delta_{ij})(\sigma_{ij}-p\delta_{ij})
=\sigma_{ij}\sigma_{ij}-2p\,\sigma_{ii}+p^2\delta_{ij}\delta_{ij}
$$

其中 $\sigma_{ii}=3p$，$\delta_{ij}\delta_{ij}=3$：

$$
s_{ij}s_{ij}=\sigma_{ij}\sigma_{ij}-6p^2+3p^2=\sigma_{ij}\sigma_{ij}-3p^2
$$

**依据**：$\delta_{ij}\delta_{ij}=\delta_{ii}=3$（三维）。

所以

$$
J_2=\tfrac12\sigma_{ij}\sigma_{ij}-\tfrac32p^2
=\tfrac12\sigma_{ij}\sigma_{ij}-\tfrac16(\sigma_{kk})^2
$$

**（c）用主应力表示。** 取主轴（$\sigma_{ij}$ 对角），$\sigma_{ij}\sigma_{ij}=\sigma_1^2+\sigma_2^2+\sigma_3^2$，
$p^2=(\sigma_1+\sigma_2+\sigma_3)^2/9$：

$$
J_2=\tfrac12(\sigma_1^2+\sigma_2^2+\sigma_3^2)-\tfrac16(\sigma_1+\sigma_2+\sigma_3)^2
$$

展开并合并：

$$
J_2=\frac{1}{6}\left[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2\right]
$$

**依据**：展开右侧得 $2(\sigma_1^2+\sigma_2^2+\sigma_3^2)-2(\sigma_1\sigma_2+\sigma_2\sigma_3+\sigma_3\sigma_1)$，
再除以 6，与上式逐项相等。✓

> 用另一种常见等价形式验算：$J_2=\frac13(\sigma_1^2+\sigma_2^2+\sigma_3^2-\sigma_1\sigma_2-\sigma_2\sigma_3-\sigma_3\sigma_1)$，两者相同。

### 2.3 第三步：为什么用 $\sqrt{3J_2}$ 作为"等效应力"

要求：找一个标量 $\sigma_{\mathrm v}(\boldsymbol{\sigma})$，使得

- 对单轴拉伸 $\sigma_{\mathrm v}=\sigma$（这样能直接与单轴屈服强度比较）；
- 对静水压力 $\sigma_{\mathrm v}=0$（静水不引起屈服）；
- 与坐标系无关（是不变量的函数）。

最自然的选择是用 $J_2$（它是"畸变强度"的最低阶非零不变量）。设 $\sigma_{\mathrm v}=c\sqrt{J_2}$。

单轴 $\sigma_1=\sigma,\sigma_2=\sigma_3=0$：

$$
J_2=\frac16\left[\sigma^2+\sigma^2+0\right]=\frac{\sigma^2}{3}
$$

要 $\sigma_{\mathrm v}=\sigma$，则 $c\sqrt{\sigma^2/3}=\sigma\Rightarrow c=\sqrt3$。

$$
\boxed{\;\sigma_{\mathrm v}=\sqrt{3J_2}\;}
$$

**依据**：标定常数 $c$ 的唯一途径就是"与单轴试验对齐"。

### 2.4 第四步：物理含义——畸变能

线弹性应变能密度可以分成两部分。先看总式：

$$
W=\tfrac12\sigma_{ij}\varepsilon_{ij}
$$

代入 $\sigma_{ij}=p\delta_{ij}+s_{ij}$ 与 $\varepsilon_{ij}=\frac{1}{3}\theta\delta_{ij}+e_{ij}$（$\theta=\varepsilon_{kk}$，$e_{ij}$ 为偏应变）：

$$
W=\tfrac12 p\,\theta+\tfrac12 s_{ij}e_{ij}
$$

**依据**：交叉项因正交而消失（与 2.1 同理）。

再由 $p=K\theta$、$s_{ij}=2Ge_{ij}$：

$$
W=\underbrace{\frac{p^2}{2K}}_{\text{体积改变}}+\underbrace{\frac{J_2}{2G}}_{\text{形状改变}}
$$

**依据**：各向同性线弹性，$K$ 体积模量、$G$ 剪切模量；$e_{ij}=s_{ij}/(2G)$。

**结论**：$\sigma_{\mathrm v}$ 相同 ⟺ $J_2$ 相同 ⟺ **畸变能相同**。
von Mises 准则"畸变能密度达到临界值即屈服"，本质就是"$\sigma_{\mathrm v}$ 达到单轴屈服强度"。

### 2.5 第五步：三个经典算例（务必背下来）

| 状态 | 主应力 (MPa) | $J_2$ | $\sigma_{\mathrm v}$ |
|---|---|---|---|
| 单轴拉伸 100 | $(100,0,0)$ | $100^2/3=3333$ | **100** |
| 纯剪 50 | $(50,0,-50)$ | $50^2=2500$ | **86.6** |
| 静水压缩 | $(-100,-100,-100)$ | $0$ | **0** |
| 等双轴拉伸 | $(100,100,0)$ | $100^2/3=3333$ | **100** |

**纯剪验算**：$\frac16[(50-0)^2+(0+50)^2+(-50-50)^2]=\frac16[2500+2500+10000]=2500$，
$\sqrt{3\times2500}=86.6$。✓

**两个反直觉的结论**：

- **纯剪 50 MPa 的"等效"程度相当于单轴 86.6 MPa**——剪切更"危险"，因为有三个主应力差同时在起作用。
- **静水 100 MPa 的等效应力是 0**——材料内部仍有 100 MPa 正应力，但不会因畸变而屈服。
  这就是为什么深海里的金属球不会被"压扁"。

> 图 4：平面应力下三种准则的屈服边界（无量纲化，$\sigma_{a,\mathrm{eq}}=1$）。
> 文件：`assets/figures/fig04-yield-criteria.png`

---

## 3. 数值检验

**检验 1｜不变量自洽**
对 $(100,40,-20)$：$J_2=\frac16[60^2+60^2+(-120)^2]=\frac16[3600+3600+14400]=3600$。
$\sigma_{\mathrm v}=\sqrt{10800}=103.9$ MPa。

用另一种写法：$\sigma_{ij}\sigma_{ij}=100^2+40^2+20^2=12000$，
$p=40$，$J_2=\frac12(12000)-\frac16(120)^2=6000-2400=3600$。✓

**检验 2｜与莫尔圆一致**
$\tau_{\max}=(\sigma_1-\sigma_3)/2=60$ MPa（见 [04](04-mohr-circle.md) §3.4）。
$\sigma_{\mathrm v}=\sqrt3\,\tau_{\max}$ 仅在 $\sigma_2=(\sigma_1+\sigma_3)/2$ 时成立；
本题 $\sigma_2=40$ 恰为中值，故 $\sigma_{\mathrm v}=\sqrt3\times60=103.9$。✓

**检验 3｜单位**
$J_2$ 的量纲是应力的平方，$\sigma_{\mathrm v}$ 才是应力。写"$J_2=3600$ MPa"是错的，
应写 $J_2=3600\ \mathrm{MPa^2}$。

---

## 4. 常见误解

**误解 1："等效应力就是最大主应力。"**

为什么错：两者数值一般不同。单轴时相等（100 vs 100），纯剪时差 1.73 倍（86.6 vs 50）。
等效应力还包含 $\sigma_2$ 的贡献，最大主应力不含。

**误解 2："静水压力下材料不受力，所以应力为零。"**

为什么错：应力不为零（三向 100 MPa），只是**偏应力为零**、等效应力为零。
"不屈服"和"没有应力"是两件事。

**误解 3："等效应力云图能看出破坏位置，所以它是真实应力。"**

为什么错：它是人为构造的标量。真实材料在缺口处可能因**最大主应力**（拉）开裂，
而等效应力的热点在别处。两个云图热点重合纯属巧合。

**误解 4："von Mises 适用于所有材料。"**

为什么错：它基于各向同性、延性、畸变能准则。
- 岩石、混凝土、铸铁：拉压强度差异大，常用 **Mohr–Coulomb / Hoek–Brown**；
- 脆性开裂：常用 **最大主应力准则（Rankine）** 或断裂力学；
- 各向异性（页岩）：需要横观各向同性准则。

---

## 5. 自测

1. 求 $(80,20,0)$ MPa 的 $p$、$J_2$、$\sigma_{\mathrm v}$。
2. 证明：$\sigma_1=\sigma_2=\sigma_3$ 时 $J_2=0$。
3. 纯剪状态下，$\sigma_{\mathrm v}/\tau_{\max}$ 是多少？单轴拉伸时又是多少？
4. 为什么"静水压力不引起屈服"对金属成立，对岩石却只是近似？（提示：摩擦、压硬、孔隙）
5. 某砂岩三轴试验：围压 20 MPa，轴压 80 MPa（压为正）。用拉正约定写出主应力，并求 $\sigma_{\mathrm v}$。

<details>
<summary>参考答案</summary>

1. $p=(80+20+0)/3=33.33$ MPa；
   $J_2=\frac16[(60)^2+(20)^2+(-80)^2]=\frac16(3600+400+6400)=1733.3\ \mathrm{MPa^2}$；
   $\sigma_{\mathrm v}=\sqrt{3\times1733.3}=72.1$ MPa。
2. 令三者均为 $\sigma$：每个差值项为 0，故 $J_2=0$。
3. 纯剪：$\sigma_{\mathrm v}=86.6$，$\tau_{\max}=50$，比值 $\sqrt3\approx1.732$。
   单轴：$\sigma_{\mathrm v}=100$，$\tau_{\max}=50$，比值 2。
4. 金属屈服由位错滑移驱动，对静水压力不敏感；
   岩石的强度来自摩擦与胶结，围压会压紧微裂隙、提高摩擦阻力，
   所以静水压力虽不直接"剪切"，却显著改变强度——需要用 $\sigma'=\sigma-\alpha p$ 与摩擦型准则。
5. 拉正：$\sigma_1=-20$（围压为压、最小），$\sigma_3=-80$（轴压、最大压）。
   $J_2=\frac16[(-20+80)^2+(-80+20)^2+0]=\frac16[3600+3600]=1200$；
   $\sigma_{\mathrm v}=\sqrt{3600}=60$ MPa。
   （也可直接用差：$\sigma_1-\sigma_3=60$，单轴压缩时 $\sigma_{\mathrm v}=|\sigma_1-\sigma_3|$。）
</details>

---

## 延伸阅读

- von Mises, R. (1913). *Mechanik der festen Körper im plastisch-deformablen Zustand*.
- Hill, R. (1950). *The Mathematical Theory of Plasticity*, Ch. 2.
- Jaeger, Cook & Zimmerman (2007). *Fundamentals of Rock Mechanics*, Ch. 4（应力与不变量）。

下一篇 → [从固体力学到地质力学：压正、地应力与有效应力](06-from-solid-to-geomechanics.md)
