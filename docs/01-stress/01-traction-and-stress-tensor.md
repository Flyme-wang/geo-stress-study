# 面力与应力张量：为什么"同一个点"还要指定截面方向

> 对应知识图谱节点：`M1-01`
> 前置：[符号与约定](../00-notation/symbols-and-conventions.md)
> 后续：[斜截面上的正应力与剪应力](02-stress-transformation.md)

---

## 0. 结论先行

1. **力**作用在物体上；**面力**（traction）作用在"过某点的某个截面"上。
   应力不是力，而是"单位面积上的面力"在一点处的极限。
2. 过同一点可以切出无穷多个截面，每个截面上的面力一般不同。
   于是"这一点的应力"不能是一个矢量，只能是**一个把方向映射成矢量的线性算子**——这就是应力张量 $\boldsymbol{\sigma}$：

$$
\mathbf{t}(\mathbf{n}) = \boldsymbol{\sigma}\,\mathbf{n}
\qquad\text{即}\qquad
t_i = \sigma_{ij} n_j
$$

3. 因此"应力需要指定截面方向"不是约定的麻烦，而是**力的矢量性 + 截面的方向性**共同决定的必然结果。
4. 在无体力偶的连续介质中，应力张量对称（6 个独立分量），并且是**二阶张量**：坐标系旋转时按 $\boldsymbol{\sigma}' = \mathbf{Q}\boldsymbol{\sigma}\mathbf{Q}^{\mathsf T}$ 变换。

适用范围：连续介质、小变形、无体力偶、准静态（惯性力可忽略）。

---

## 1. 前提与符号

| 假设 | 内容 | 不满足会怎样 |
|---|---|---|
| H1 连续介质 | 物质点连续分布，可定义"面积趋于零"的极限 | 在裂纹面、界面、孔隙尺度上失效，需用接触/断裂力学 |
| H2 无体力偶 | 物质微元上不作用分布力偶 | 应力张量不对称，出现 Cosserat 介质 |
| H3 准静态 | 惯性项可忽略（或已并入体力） | 需写 $\rho\ddot{u}$ 项，结论形式不变但不再是静力平衡 |
| H4 小变形 | 用变形前构型计算面积与法向 | 大变形需用 Cauchy 应力与 Kirchhoff 应力区分 |

符号：

- $\Delta S$：截面面积；$\Delta \mathbf{F}$：该截面上的内力合力。
- $\mathbf{n}$：截面单位法向量，$n_i n_i = 1$。
- $\mathbf{t}(\mathbf{n})$：法向为 $\mathbf{n}$ 的截面上的面力矢量，单位 Pa。
- $\sigma_{ij}$：应力张量分量，$i$ 为**作用面法向**，$j$ 为**力的方向**。

---

## 2. 推导

### 2.1 第一步：定义面力（traction）

取过点 $P$、法向为 $\mathbf{n}$ 的一个小面元，面积为 $\Delta S$，面上的内力合力为 $\Delta \mathbf{F}$。定义

$$
\mathbf{t}(\mathbf{n}) \;=\; \lim_{\Delta S \to 0}\frac{\Delta \mathbf{F}}{\Delta S}
$$

**依据**：这是"应力"的唯一严格定义方式。极限存在依赖 H1。

**注意两点**：

1. 分母是**面积**，不是体积，也不是"物体总横截面积"。
2. 面力是矢量，有三个分量；普通意义的"应力"往往指它的某个分量（正应力或剪应力）。

> 单位：$1\ \mathrm{Pa} = 1\ \mathrm{N/m^2}$；$1\ \mathrm{MPa} = 1\ \mathrm{N/mm^2}$。

### 2.2 第二步：为什么必须给出 $\mathbf{n}$

沿不同方向切同一个点，切出来的面承受的内力一般不同——因为两侧材料"拉住的面积分布"不同。
所以 $\mathbf{t}$ 是 $\mathbf{n}$ 的函数：$\mathbf{t} = \mathbf{t}(\mathbf{n})$。

若只写"$P$ 点应力是 100 MPa"，信息是不完整的：它没说这个 100 MPa 是哪个面上的、沿哪个方向。

### 2.3 第三步：反作用——$\mathbf{t}(-\mathbf{n}) = -\mathbf{t}(\mathbf{n})$

取一块极薄的"平板"微元，两个面法向分别为 $\mathbf{n}$ 与 $-\mathbf{n}$，面积均为 $\Delta S$，厚度 $h$。

对这块板列力平衡：

$$
\mathbf{t}(\mathbf{n})\Delta S + \mathbf{t}(-\mathbf{n})\Delta S + \rho\mathbf{b}\,\Delta S\,h = \mathbf{0}
$$

**依据**：牛顿第二定律（准静态，合力为零）。

两边除以 $\Delta S$，令 $h \to 0$：

$$
\mathbf{t}(-\mathbf{n}) = -\,\mathbf{t}(\mathbf{n})
$$

这一步只用到了**厚度趋于零时体积力是高阶小量**，是后面所有推导的关键技巧。

### 2.4 第四步：四面体论证 → Cauchy 应力定理

在点 $P$ 附近取一个四面体：

- 三个面分别平行于三个坐标面，外法向为 $-\mathbf{e}_1, -\mathbf{e}_2, -\mathbf{e}_3$，面积记为 $\Delta S_1,\Delta S_2,\Delta S_3$；
- 第四个斜面外法向为 $\mathbf{n}$，面积 $\Delta S$；
- 四面体高 $h$（$P$ 到斜面的距离）。

**几何关系（投影公式）**：斜面的面积投影到坐标面上，得到三个直角面：

$$
\Delta S_i = n_i\,\Delta S
$$

**依据**：面积投影定理——一个平面面积在法向为 $\mathbf{e}_i$ 的平面上的投影，等于原面积乘以两法向夹角余弦，而 $\cos\angle(\mathbf{n},\mathbf{e}_i)=n_i$。
可自行验证：$n_i$ 满足 $n_1^2+n_2^2+n_3^2=1$，三个投影面积平方和等于 $\Delta S^2$，与几何一致。

**力平衡**：

$$
\mathbf{t}(\mathbf{n})\,\Delta S
+\sum_{i=1}^{3}\mathbf{t}(-\mathbf{e}_i)\,\Delta S_i
+\rho\mathbf{b}\,\Delta V = \mathbf{0},
\qquad \Delta V = \tfrac13 h\,\Delta S
$$

**依据**：四面体整体静力平衡（H3）。

用 2.3 步的结论 $\mathbf{t}(-\mathbf{e}_i) = -\mathbf{t}(\mathbf{e}_i)$ 替换，并代入 $\Delta S_i = n_i\Delta S$：

$$
\mathbf{t}(\mathbf{n})\,\Delta S
-\sum_{i=1}^{3}\mathbf{t}(\mathbf{e}_i)\,n_i\,\Delta S
+\rho\mathbf{b}\,\tfrac13 h\,\Delta S = \mathbf{0}
$$

两边除以 $\Delta S$：

$$
\mathbf{t}(\mathbf{n}) = \sum_{i=1}^{3} n_i\,\mathbf{t}(\mathbf{e}_i) \;-\; \tfrac13\rho h\,\mathbf{b}
$$

**令 $h\to 0$**（把四面体缩到点 $P$），体积力项消失：

$$
\boxed{\;\mathbf{t}(\mathbf{n}) = n_1\mathbf{t}(\mathbf{e}_1)+n_2\mathbf{t}(\mathbf{e}_2)+n_3\mathbf{t}(\mathbf{e}_3)\;}
$$

这一步的结论极其重要：**任意方向截面上的面力，都可以由三个坐标面上的面力线性组合得到，权重就是法向分量 $n_i$。**
换句话说，$\mathbf{t}$ 对 $\mathbf{n}$ 是**线性**的。这不是假设，是由平衡推出的定理（Cauchy 应力定理）。

### 2.5 第五步：引入应力张量

把三个坐标面上的面力按坐标轴分解，定义 9 个数：

$$
\mathbf{t}(\mathbf{e}_i) = \sigma_{i1}\mathbf{e}_1 + \sigma_{i2}\mathbf{e}_2 + \sigma_{i3}\mathbf{e}_3
\qquad\Longleftrightarrow\qquad
t_j(\mathbf{e}_i) = \sigma_{ij}
$$

即 $\sigma_{ij}$ = "法向为 $\mathbf{e}_i$ 的面上，沿 $\mathbf{e}_j$ 方向的面力分量"。

代入 2.4 的结果：

$$
t_j(\mathbf{n}) = \sum_{i=1}^{3}\sigma_{ij}\,n_i
\qquad\Longleftrightarrow\qquad
\boxed{\;t_i = \sigma_{ij}n_j\ \ \text{或}\ \ \mathbf{t} = \boldsymbol{\sigma}^{\mathsf T}\mathbf{n}\;}
$$

写成矩阵：

$$
\begin{bmatrix}t_1\\t_2\\t_3\end{bmatrix}
=
\begin{bmatrix}
\sigma_{11} & \sigma_{12} & \sigma_{13}\\
\sigma_{21} & \sigma_{22} & \sigma_{23}\\
\sigma_{31} & \sigma_{32} & \sigma_{33}
\end{bmatrix}^{\mathsf T}
\begin{bmatrix}n_1\\n_2\\n_3\end{bmatrix}
$$

**9 个分量的含义**（以图 1 的微元为例）：

| 分量 | 所在面 | 方向 | 类型 |
|---|---|---|---|
| $\sigma_{11},\sigma_{22},\sigma_{33}$ | 垂直各坐标轴的面 | 沿该面法向 | 正应力 |
| $\sigma_{12},\sigma_{13},\sigma_{21},\sigma_{23},\sigma_{31},\sigma_{32}$ | 同上 | 在面内 | 剪应力 |

> 图 1：三维微元的三个可见面上，绿色箭头为正应力分量，红色箭头为剪应力分量。
> 注意每个剪应力都有下标：第一个指标是"在哪个面上"，第二个是"指向哪个方向"。
> 文件：`assets/figures/fig01-stress-cube.png`

### 2.6 第六步：为什么 9 个变 6 个——对称性

对边长 $a$ 的立方体微元，关于 $x_3$ 轴取力矩平衡。

作用在 $\pm x_1$ 面上的剪应力 $\sigma_{12}$ 形成一对力偶，力臂 $a$，力的大小 $\sigma_{12}a^2$；
作用在 $\pm x_2$ 面上的 $\sigma_{21}$ 同理。取矩（逆时针为正）：

$$
M_3 = \underbrace{\sigma_{12}\,a^2\cdot a}_{\text{一对力偶}} \;-\; \underbrace{\sigma_{21}\,a^2\cdot a}_{\text{反向一对}} \;+\; \underbrace{\rho b\,a^3\cdot \tfrac{a}{2}}_{\text{体积力}} \;=\; 0
$$

**依据**：角动量守恒（准静态，角加速度为零；H2 保证没有分布力偶）。

两边除以 $a^3$：

$$
\sigma_{12}-\sigma_{21} + \tfrac12\rho b\,a = 0
$$

令 $a\to 0$，体积力项消失：

$$
\sigma_{12}=\sigma_{21}
$$

同理 $\sigma_{13}=\sigma_{31}$、$\sigma_{23}=\sigma_{32}$。

$$
\boxed{\;\sigma_{ij}=\sigma_{ji}\;}
$$

于是独立分量从 9 个降为 **6 个**：

$$
(\sigma_{11},\ \sigma_{22},\ \sigma_{33},\ \sigma_{12},\ \sigma_{13},\ \sigma_{23})
$$

> 这一条是"能量守恒/角动量守恒"的结果，不是定义。若介质能承受分布力偶（如某些颗粒材料模型），对称性不成立。

### 2.7 第七步：为什么是"张量"——坐标变换

设新坐标系 $\{\mathbf{e}'_i\}$ 由旧坐标系经正交旋转 $\mathbf{Q}$ 得到，$Q_{ij} = \mathbf{e}'_i\cdot\mathbf{e}_j$，满足 $\mathbf{Q}\mathbf{Q}^{\mathsf T}=\mathbf{I}$。

同一个面力矢量在两套坐标系下的分量关系为 $t'_i = Q_{ij}t_j$；同一个法向为 $n'_i = Q_{ij}n_j$。

由 $t'_i = \sigma'_{ij}n'_j$ 与 $t_i=\sigma_{ij}n_j$：

$$
\sigma'_{ij}n'_j = Q_{ij}t_j = Q_{ij}\sigma_{jk}n_k = Q_{ij}\sigma_{jk}Q_{lk}n'_l
$$

由于此式对任意 $\mathbf{n}'$ 成立，且 $Q_{lk} = (Q^{\mathsf T})_{kl}$：

$$
\boxed{\;\sigma'_{ij} = Q_{ik}\,Q_{jl}\,\sigma_{kl}
\qquad\Longleftrightarrow\qquad
\boldsymbol{\sigma}' = \mathbf{Q}\,\boldsymbol{\sigma}\,\mathbf{Q}^{\mathsf T}\;}
$$

**依据**：正交变换的定义与"矢量分量随坐标系旋转"的规则。

"满足这种变换律的 9 个量"就是**二阶张量**的定义。它跟"随便找 9 个数填进矩阵"的区别在于：换坐标系时，这 9 个数**必须**按上式变。

二维特例（绕 $x_3$ 转 $\theta$）：

$$
\mathbf{Q} =
\begin{bmatrix}
\cos\theta & \sin\theta\\
-\sin\theta & \cos\theta
\end{bmatrix}
$$

代入即得下一节的斜截面公式（见 [02](02-stress-transformation.md)）。

---

## 3. 数值检验

**检验 1：单位与量级**
$100\ \mathrm{MPa}$ 的拉应力作用在 $1\ \mathrm{cm^2}$ 截面上，内力为
$100\times10^6 \times 10^{-4} = 10^4\ \mathrm{N} \approx 1\ \mathrm{ton}$ 力。量级合理。

**检验 2：单轴状态的面力**

取 $\sigma_{11}=100$ MPa，其余为 0。

- 法向 $\mathbf{n}=(1,0,0)$（垂直于棒轴的截面）：
  $t_1=\sigma_{11}n_1=100$ MPa，$t_2=t_3=0$。面力全为法向，正应力 100 MPa，剪应力 0。✓
- 法向 $\mathbf{n}=(0,1,0)$（平行于棒轴的截面）：
  $t_i=\sigma_{i2}n_2=0$。**自由表面，无面力**。✓
- 法向 $\mathbf{n}=(1,1,0)/\sqrt2$（45° 斜截面）：
  $t_1 = 100\times\tfrac{1}{\sqrt2} = 70.71$ MPa，$t_2=t_3=0$。
  面力大小 70.71 MPa，方向仍沿 $x_1$。
  正应力 $\sigma_n = n_i t_i = 70.71/\sqrt2 = 50$ MPa；
  剪应力 $\tau=\sqrt{70.71^2-50^2}=50$ MPa。✓

**检验 3：三角恒等式自洽**
$\sigma_n^2+\tau^2 = |\mathbf{t}|^2$。
$50^2+50^2 = 5000 = 70.71^2$。✓

**检验 4：对称性量纲**
力矩平衡式中 $\sigma a^3$ 与 $\rho b a^4$ 同量纲（N·m）：
Pa·m³ = N/m²·m³ = N·m；kg/m³·m/s²·m⁴ = N·m。✓

---

## 4. 常见误解

**误解 1："某点的应力 = 外力 / 截面积。"**

为什么错：$\mathbf{F}/A$ 只有在**均匀受力且截面固定**时才等于该截面上的平均面力；
一点处的应力是 $\Delta S\to0$ 的极限，与外载荷总量不是一一对应。同一根棒，同一点，不同方向的截面上面力不同。

正确说法：外力决定**边界条件**，应力是内部逐点分布的场，靠平衡方程与边界条件求解。

**误解 2："应力是矢量，所以有三个分量。"**

为什么错：把面力（矢量）与应力（张量）混为一谈。
面力 $\mathbf{t}$ 是矢量，但它依赖 $\mathbf{n}$；要描述"所有可能截面的面力"，需要 9 个数（对称后 6 个）。

正确说法：**应力是二阶张量；面力是矢量，是应力张量与法向的乘积。**

**误解 3："没有外力，就没有应力。"**

为什么错：自重、温度、孔隙压力、相邻材料约束都会产生应力。
地壳中一块"什么都没压"的岩石，在 2000 m 深处仍有约 49 MPa 的垂向压应力。

**误解 4："$\sigma_{12}$ 和 $\sigma_{21}$ 是两个独立的剪应力。"**

为什么错：角动量守恒要求它们相等。它们描述的是同一个"剪切对"的两个侧面。

---

## 5. 自测

1. 某点应力张量（MPa）为
   $\begin{bmatrix}100&0&0\\0&0&0\\0&0&0\end{bmatrix}$。
   求法向 $\mathbf{n}=(0,0,1)$ 的截面上的面力矢量。
2. 接上题，求 $\mathbf{n}=(1,0,0)$ 面上沿 $x_1$ 的面力分量。
3. 说明为什么"过一点有无穷多个截面"不意味着"无穷多个独立应力分量"。
4. 若某介质能承受分布体力偶，第 2.6 步的推导会失效在哪一步？结论如何改变？
5. 一个学生说："这块岩石没受力，因为周围没有东西压它。"请用本页内容纠正，并给出一个数量级估计。

<details>
<summary>参考答案</summary>

1. $t_i=\sigma_{ij}n_j$，只有 $n_3=1$：$t_1=0,t_2=0,t_3=0$。该面是自由面，无面力。
2. $n_1=1$，$t_1=\sigma_{11}=100$ MPa，其余为 0。
3. 每个截面上的面力都由同一组 9 个（对称后 6 个）张量分量**线性组合**得到：
   $t_i=\sigma_{ij}n_j$。无穷多个方向对应无穷多组 $(n_1,n_2,n_3)$，但底层的独立参数只有 6 个。
4. 失效在"力矩平衡中体积力项趋于零"这一推理链上的前提：若存在分布力偶 $m$，
   平衡式变为 $\sigma_{12}-\sigma_{21}+m=0$，不再得到对称。此时张量非对称，需 9 个分量。
5. 应力来自**自重与围压**。取 $\rho=2500\ \mathrm{kg/m^3}$，$z=1\ \mathrm{km}$：
   $\sigma_v\approx 2500\times9.81\times1000 \approx 24.5$ MPa。深处岩石处于三向受压状态，
   只是"没有额外载荷"，不是"没有应力"。
</details>

---

## 延伸阅读

- Cauchy, A.-L. (1827). *Exercices de mathématiques* —— 应力概念的最初形式。
- Timoshenko, S. P. (1953). *History of Strength of Materials*，第 1 章。
- 任何一本连续介质力学教材的 "Cauchy stress theorem" 一节（如 Gurtin, *An Introduction to Continuum Mechanics*, §14）。

下一篇 → [斜截面上的正应力与剪应力](02-stress-transformation.md)
