# 主应力：切应力为什么会消失，主应力又是什么

> 对应知识图谱节点：`M1-03`
> 前置：[斜截面变换](02-stress-transformation.md)
> 后续：[莫尔圆](04-mohr-circle.md)、[偏应力与等效应力](05-deviator-and-von-mises.md)

---

## 0. 结论先行

1. 主应力是**应力张量的特征值**（eigenvalue）。求主应力 = 解一个 $3\times3$ 对称矩阵的特征值问题，不是"另外加一股力"。
2. 主方向是特征向量方向；三个主方向**互相垂直**（对称张量的谱定理）。
3. 在任意截面上，正应力 $\sigma_n$ 的最大值就是 $\sigma_1$，最小值就是 $\sigma_3$，它们出现在主截面上。这就是"最大主应力"这个说法的严格含义。
4. 三个主应力的组合

$$
I_1=\sigma_1+\sigma_2+\sigma_3,\quad
I_2=\sigma_1\sigma_2+\sigma_2\sigma_3+\sigma_3\sigma_1,\quad
I_3=\sigma_1\sigma_2\sigma_3
$$

与坐标系无关，称为**应力不变量**。

适用范围：同前一篇（连续介质、小变形、无体力偶）。

---

## 1. 前提与符号

- $\boldsymbol{\sigma}$ 为实对称二阶张量（对称性来自角动量守恒，见 [01](01-traction-and-stress-tensor.md) §2.6）。
- $\mathbf{n}$ 为单位法向， $|\mathbf{n}|=1$。
- 主应力 $\lambda$，主方向 $\mathbf{n}$，满足 $\mathbf{t}=\lambda\mathbf{n}$。
- 排序约定： $\sigma_1\ge\sigma_2\ge\sigma_3$（**代数**大小，拉正）。

---

## 2. 推导

### 2.1 第一步：主方向的定义

"某截面上没有剪应力"等价于"面力矢量完全沿法向"：

$$
\mathbf{t}(\mathbf{n}) = \lambda\ \mathbf{n}
\qquad\text{其中}\ \lambda=\sigma_n
$$

**依据**：若 $\tau=0$，面力只有法向分量，其大小即正应力，方向即 $\pm\mathbf{n}$。

代入 $\mathbf{t}=\boldsymbol{\sigma}\mathbf{n}$：

$$
\boxed{\boldsymbol{\sigma}\mathbf{n}=\lambda\mathbf{n}}
\qquad\Longleftrightarrow\qquad
(\sigma_{ij}-\lambda\delta_{ij})n_j=0
$$

这就是标准的**特征值问题**。

### 2.2 第二步：非零解的条件 → 特征方程

$\mathbf{n}=\mathbf{0}$ 无意义，所以矩阵 $(\boldsymbol{\sigma}-\lambda\mathbf{I})$ 必须奇异：

$$
\det(\boldsymbol{\sigma}-\lambda\mathbf{I})=0
$$

**依据**：线性代数基本定理——齐次方程有非零解当且仅当系数矩阵行列式为零。

三维展开（这一步只做代数，不跳步）：

$$
\det
\begin{bmatrix}
\sigma_{11}-\lambda & \sigma_{12} & \sigma_{13}\\
\sigma_{12} & \sigma_{22}-\lambda & \sigma_{23}\\
\sigma_{13} & \sigma_{23} & \sigma_{33}-\lambda
\end{bmatrix}=0
$$

按第一行展开并合并同类项，得到**特征方程**：

$$
\boxed{\lambda^3-I_1\lambda^2+I_2\lambda-I_3=0}
$$

其中

$$
I_1=\sigma_{11}+\sigma_{22}+\sigma_{33}=\sigma_{kk}
$$

$$
I_2=\sigma_{11}\sigma_{22}+\sigma_{22}\sigma_{33}+\sigma_{33}\sigma_{11}
-\sigma_{12}^2-\sigma_{23}^2-\sigma_{31}^2
=\tfrac12\left[(\sigma_{kk})^2-\sigma_{ij}\sigma_{ij}\right]
$$

$$
I_3=\det(\sigma_{ij})
$$

**依据**：行列式展开 + 对称性 $\sigma_{ij}=\sigma_{ji}$（否则 $I_2$ 中还会出现 $\sigma_{12}\sigma_{21}$ 的交叉项，形式不同）。

### 2.3 第三步：为什么三个根一定存在且为实数

这里不做完整证明，但要说清逻辑链，避免"凭感觉"：

**（a） $\sigma_n(\mathbf{n})$ 必有最大值与最小值。**
$\sigma_n(n_i n_i=1)=n_i\sigma_{ij}n_j$ 是单位球面（有界闭集）上的连续函数，
由**极值定理**必取到最大、最小值。

**（b）极值点处的法向就是主方向。**
用拉格朗日乘子处理约束 $n_in_i=1$：

$$
\frac{\partial}{\partial n_k}\Big[n_i\sigma_{ij}n_j-\lambda\ (n_in_i-1)\Big]=0
\quad\Longrightarrow\quad
2\sigma_{kj}n_j-2\lambda n_k=0
$$

即 $\boldsymbol{\sigma}\mathbf{n}=\lambda\mathbf{n}$。**依据**：约束极值的必要条件。

两边左乘 $n_k$ 并求和： $\lambda=n_k\sigma_{kj}n_j=\sigma_n(\mathbf{n})$，是**实数**。
**依据**：实数的二次型是实数。于是至少有一个实特征值，且它就是正应力的极值。

**（c）降维得到另外两个。**
在垂直于该主方向的平面内重复同样的论证（二维问题），得到一个二次方程，判别式非负
（对称矩阵 $2\times2$ 的特征值是 $\frac{\sigma_{11}+\sigma_{22}}{2}\pm\sqrt{(\frac{\sigma_{11}-\sigma_{22}}{2})^2+\sigma_{12}^2}$，根号内非负）。
第三维再降一次。于是三个实根都存在。

> 严谨的完整证明（实对称矩阵可正交对角化）属于线性代数，结论可直接引用：
> 实对称矩阵的特征值全为实数，且存在一组标准正交特征向量。

### 2.4 第四步：三个主方向互相垂直

设 $\boldsymbol{\sigma}\mathbf{n}^{(1)}=\lambda_1\mathbf{n}^{(1)}$， $\boldsymbol{\sigma}\mathbf{n}^{(2)}=\lambda_2\mathbf{n}^{(2)}$， $\lambda_1\ne\lambda_2$。

$$
\lambda_1\ \mathbf{n}^{(1)}\cdot\mathbf{n}^{(2)}
=(\boldsymbol{\sigma}\mathbf{n}^{(1)})\cdot\mathbf{n}^{(2)}
=\mathbf{n}^{(1)}\cdot(\boldsymbol{\sigma}^{\mathsf T}\mathbf{n}^{(2)})
=\mathbf{n}^{(1)}\cdot(\boldsymbol{\sigma}\mathbf{n}^{(2)})
=\lambda_2\ \mathbf{n}^{(1)}\cdot\mathbf{n}^{(2)}
$$

**依据**：点积的对称性 + 张量对称性 $\boldsymbol{\sigma}^{\mathsf T}=\boldsymbol{\sigma}$。

移项：

$$
(\lambda_1-\lambda_2)\ \mathbf{n}^{(1)}\cdot\mathbf{n}^{(2)}=0
\quad\Longrightarrow\quad
\mathbf{n}^{(1)}\cdot\mathbf{n}^{(2)}=0
$$

**结论**：不同主应力对应的主方向正交。

> 若两个特征值相等（例如静水状态三个都相等，或轴对称状态两个相等），主方向不唯一，
> 可在该子空间内任选一组正交方向。此时"主方向"是一个**面**或**体**，不是一条线。

### 2.5 第五步：主应力就是正应力的极值

由 2.3(b)，极值处的 $\lambda=\sigma_n$。
因为是三次多项式，三个根 $\lambda_1\ge\lambda_2\ge\lambda_3$ 分别是 $\sigma_n$ 在球面上的**全局最大、鞍点、全局最小**。

$$
\min_{\mathbf{n}}\sigma_n=\sigma_3,\qquad
\max_{\mathbf{n}}\sigma_n=\sigma_1
$$

这就是"最大主应力"的严格含义：**任意截面上正应力的最大值**，而不是"最大的那个外力"。

> 注意： $\sigma_2$ 是鞍点，对应的是"次大值方向"，在三维中它既不是最大也不是最小。

### 2.6 第六步：不变量为什么"不变"

由张量变换律 $\boldsymbol{\sigma}'=\mathbf{Q}\boldsymbol{\sigma}\mathbf{Q}^{\mathsf T}$，
特征多项式

$$
\det(\boldsymbol{\sigma}'-\lambda\mathbf{I})
=\det(\mathbf{Q}(\boldsymbol{\sigma}-\lambda\mathbf{I})\mathbf{Q}^{\mathsf T})
=\det(\mathbf{Q})\det(\boldsymbol{\sigma}-\lambda\mathbf{I})\det(\mathbf{Q}^{\mathsf T})
=\det(\boldsymbol{\sigma}-\lambda\mathbf{I})
$$

**依据**：行列式的乘法性质 + $\det\mathbf{Q}=\pm1$（旋转取 $+1$）。

多项式系数不变，所以 $I_1,I_2,I_3$ 与坐标系无关。它们可以写成主应力的组合：

$$
I_1=\sigma_1+\sigma_2+\sigma_3
$$

$$
I_2=\sigma_1\sigma_2+\sigma_2\sigma_3+\sigma_3\sigma_1
$$

$$
I_3=\sigma_1\sigma_2\sigma_3
$$

---

## 3. 数值检验

### 3.1 单向拉伸

$\sigma_{11}=100$，其余为 0。特征方程
$\lambda^3-100\lambda^2=0\Rightarrow \lambda=100,0,0$。

主应力 $(100,0,0)$ MPa。主方向： $x_1$（任意垂直于它的方向也是主方向，因为重根）。

**交叉验证**：用 [02](02-stress-transformation.md) 的公式， $\theta=0$ 时 $\tau=0$、 $\sigma_n=100$；
$\theta=90^\circ$ 时 $\tau=0$、 $\sigma_n=0$。✓

### 3.2 纯剪 50 MPa

$\tau_{12}=50$，正应力分量为 0。特征方程：

$$
\lambda^3+I_2\lambda=0,\qquad I_2=-50^2=-2500
$$

$$
\lambda(\lambda^2-2500)=0\quad\Longrightarrow\quad\lambda=50,\ 0,\ -50
$$

主应力 $(50,0,-50)$ MPa。主方向为 $x_1$ 与 $x_2$ 之间的 $45^\circ$ 方向。

**交叉验证**： $\theta=45^\circ$ 时
$\sigma_n=0+\frac{0-0}{2}\cos90^\circ+50\sin90^\circ=50$， $\tau=50\cos90^\circ=0$。✓

### 3.3 一般二维状态

$\sigma_x=80$， $\sigma_y=20$， $\tau_{xy}=40$ MPa。

二维特征值公式：

$$
\lambda_{1,2}=\frac{\sigma_x+\sigma_y}{2}\pm\sqrt{\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2}
=50\pm\sqrt{30^2+40^2}=50\pm50
$$

$$
\sigma_1=100,\quad \sigma_2=0,\quad \sigma_3=0\ (\text{平面应力})
$$

主方向：

$$
\tan 2\theta_p=\frac{2\tau_{xy}}{\sigma_x-\sigma_y}=\frac{80}{60}
\quad\Longrightarrow\quad
2\theta_p=53.13^\circ,\quad \theta_p=26.57^\circ
$$

**检验**：把 $\theta=26.57^\circ$ 代入 [02](02-stress-transformation.md) 的正应力公式：

$$
\sigma_n=50+30\cos53.13^\circ+40\sin53.13^\circ
=50+30(0.6)+40(0.8)=50+18+32=100\ \checkmark
$$

$$
\tau=-30\sin53.13^\circ+40\cos53.13^\circ=-30(0.8)+40(0.6)=0\ \checkmark
$$

### 3.4 静水应力

所有 $\sigma_n=100$，剪应力为 0。特征方程三重根 $\lambda=100$。
**任意方向都是主方向**，任意截面都是主截面。这也解释了为什么静水压力不引起塑性屈服（见下一篇）。

---

## 4. 常见误解

**误解 1："主应力是作用在物体上的另一个力。"**

为什么错：主应力没有对应的外力。它只是"把同一张量换到一组特殊坐标轴上的三个对角分量"。
外力、边界条件完全不变，只是观察的截面方向变了。

正确说法：主应力 = 特征值 = 正应力极值；主方向 = 特征向量方向。

**误解 2："最大主应力就是最大剪应力所在的截面上的正应力。"**

为什么错：最大剪应力出现在与 $\sigma_1$、 $\sigma_3$ 都成 $45^\circ$ 的截面上，那里正应力是 $(\sigma_1+\sigma_3)/2$，
**不是** $\sigma_1$。两个极值发生在不同截面上。

**误解 3："三个主应力排序时，压应力大的排前面。"**

为什么错：本库用拉正， $\sigma_1\ge\sigma_2\ge\sigma_3$ 是**代数**比较。
在压正约定下顺序整体反转，见 [符号约定](../00-notation/symbols-and-conventions.md) §5。
混用会导致莫尔圆画反、莫尔–库仑公式符号错。

**误解 4："三维里最大正应力一定作用在最大主应力方向，所以别的方向不重要。"**

为什么错：正应力最大 ≠ 剪应力最大 ≠ 能量密度最大 ≠ 最先破坏。
岩石抗压远强于抗拉，破坏常由**最小主应力（最大压）**或**剪应力**控制，取决于准则。

---

## 5. 自测

1. 求 $\sigma_{11}=60,\sigma_{22}=60,\sigma_{33}=-20$ MPa 的三个主应力与 $I_1$。
2. 某点 $\sigma_x=0,\sigma_y=0,\tau_{xy}=-50$ MPa。主应力是多少？主方向角度？
3. 证明：若 $\sigma_1=\sigma_2\ne\sigma_3$，主方向不唯一，描述其几何形状。
4. 用 2.5 的结论解释：为什么"最大主应力云图"能回答"哪里受拉最厉害"？
5. 一个同学说"主应力为 100、0、-100，所以最大剪应力是 100"。错在哪里？

<details>
<summary>参考答案</summary>

1. 对角矩阵，主应力即对角元： $(60,60,-20)$； $I_1=60+60-20=100$ MPa。
2. $\lambda^2-2500=0\Rightarrow \lambda=\pm50$（纯剪， $\sigma_3=0$ 因平面应力）→ $(50,0,-50)$。
   $\tan2\theta_p = 2(-50)/0\to\infty$， $2\theta_p=-90^\circ$， $\theta_p=-45^\circ$。
3. $\sigma_1=\sigma_2$ 时，垂直于 $\sigma_3$ 主方向的整个平面内任意方向都满足特征方程（退化）。
   主方向集合 = 一个平面内的所有方向； $\sigma_3$ 方向仍唯一。
4. 因为 $\sigma_1=\max_{\mathbf{n}}\sigma_n$：它就是在所有可能的截面上正应力的上限，
   所以云图上的 $\sigma_1$ 分布直接给出"最容易被拉裂"的位置与方向信息。
5. 最大剪应力是 $(\sigma_1-\sigma_3)/2=(100-(-100))/2=100$ MPa，数值对，
   但他把 $\sigma_1-\sigma_3$ 说成了" $\sigma_1$ 本身"。若状态是 $(100,0,0)$，最大剪应力是 50 而不是 100。
</details>

---

下一篇 → [莫尔圆：把变换公式画成一个圆](04-mohr-circle.md)
