# 符号、单位与约定

> 这一页是全库的"词典"。后面的推导只在这里定义一次符号，其余地方直接使用。
> 如果你只读一页，就读这一页——绝大多数"看不懂应力"的问题，根源是符号约定混用。

---

## 1. 三条基本约定

| 约定 | 本库采用 | 常见冲突 |
|---|---|---|
| 正负号 | **拉为正**（固体力学） | 岩石力学/地质学常用**压为正** |
| 主应力排序 | $\sigma_1 \ge \sigma_2 \ge \sigma_3$（代数值） | 有的教材按 $\sigma_1$ 为最大压应力排序 |
| 张量指标顺序 | $\sigma_{ij}$： $i$ = 作用面法向， $j$ = 力方向 | 少数教材用相反顺序（不影响对称张量的数值） |

换算规则见第 5 节。**在任何一页里，先确认这一页用的是哪套约定，再读公式。**

---

## 2. 应力相关符号

| 符号 | 名称 | 定义 / 说明 |
|---|---|---|
| $\mathbf{t}(\mathbf{n})$ | 面力矢量（traction） | 单位面积上的内力，单位 Pa |
| $n_i$ | 单位法向 | $n_i n_i = 1$ |
| $\sigma_{ij}$ | 应力张量 | $\mathbf{t} = \boldsymbol{\sigma}\cdot\mathbf{n}$，即 $t_i = \sigma_{ij}n_j$ |
| $\sigma_n$ | 正应力 | 面力沿法向的投影： $\sigma_n = n_i \sigma_{ij} n_j$ |
| $\tau$ | 剪应力 | 面力在截面内的分量： $\tau^2 = \vert\mathbf{t}\vert^2 - \sigma_n^2$ |
| $\sigma_1,\sigma_2,\sigma_3$ | 主应力 | $\sigma_{ij}$ 的三个特征值 |
| $I_1,I_2,I_3$ | 应力不变量 | 见主应力一页，与坐标系无关 |
| $p$ | 平均正应力 | $p = I_1/3 = (\sigma_1+\sigma_2+\sigma_3)/3$ |
| $s_{ij}$ | 偏应力张量 | $s_{ij} = \sigma_{ij} - p\ \delta_{ij}$ |
| $J_2$ | 偏应力第二不变量 | $J_2 = \tfrac12 s_{ij}s_{ij}$ |
| $\sigma_{\mathrm{v}}$ | von Mises 等效应力 | $\sigma_{\mathrm{v}} = \sqrt{3J_2}$ |
| $\theta$ | 截面倾角 | 从 $x_1$ 轴逆时针量到截面**法线** |
| $R$ | 莫尔圆半径 | $R = (\sigma_1-\sigma_3)/2$ |

约定 $\theta$ 的方向很重要：本库所有二维公式都规定 **$\theta$ 是法线相对 $x_1$ 轴的夹角**，
所以 $\theta=0$ 对应"法线沿 $x_1$"的截面，也就是垂直于 $x_1$ 轴的那个面。

---

## 3. 变形与强度符号

| 符号 | 名称 | 说明 |
|---|---|---|
| $\varepsilon_{ij}$ | 小应变张量 | $\varepsilon_{ij} = \tfrac12(\partial_j u_i + \partial_i u_j)$ |
| $E,\ \nu$ | 弹性模量、泊松比 | 各向同性线弹性 |
| $G$ | 剪切模量 | $G = E/[2(1+\nu)]$ |
| $K$ | 体积模量 | $K = E/[3(1-2\nu)]$ |
| $c,\ \varphi$ | 黏聚力、内摩擦角 | Mohr–Coulomb 准则 |
| $\sigma_t$ | 抗拉强度 | 拉为正 |
| $K_I$ | I 型应力强度因子 | 单位 $\mathrm{Pa}\sqrt{\mathrm{m}}$ |
| $G_c$ | 临界能量释放率 | 单位 $\mathrm{J/m^2}$ |

---

## 4. 单位换算

| 量 | 换算 |
|---|---|
| 1 MPa | $10^6$ Pa = 1 N/mm² |
| 1 GPa | $10^3$ MPa |
| 1 bar | 0.1 MPa |
| 1 psi | 6.8948 kPa |
| 1 kgf/cm² | 0.09807 MPa |
| 密度 × 重力 × 深度 | $\sigma_v \approx \rho g z$，取 $\rho=2500\ \mathrm{kg/m^3}$ 时约 **0.0245 MPa/m** |

**量纲自检法**：应力乘面积得力（N），应力乘体积得能量（J）。
推导中如果出现"应力 + 力"这样的组合，一定错了。

---

## 5. 拉正 ↔ 压正 换算

设 $\sigma^{\mathrm{T}}$ 为"拉为正"的应力张量， $\sigma^{\mathrm{C}}$ 为"压为正"的应力张量：

$$
\sigma^{\mathrm{C}}_{ij} = -\ \sigma^{\mathrm{T}}_{ij}
\qquad\Longleftrightarrow\qquad
\sigma^{\mathrm{T}}_{ij} = -\ \sigma^{\mathrm{C}}_{ij}
$$

推论（**极易出错，务必记住**）：

1. 主应力**数值**全部反号，**排序顺序反转**：
   拉正下的 $(\sigma_1,\sigma_2,\sigma_3) = (10,0,-30)$ MPa，
   在压正下变成 $(30,0,-10)$ MPa。同一个状态，写成两组完全不同的数。
2. 剪应力在正负号统一反号时**不变号**（因为它是两个方向的乘积）。
3. 莫尔圆**整体左右翻转**，圆的形状不变。
4. 有效应力公式写法不同：拉正 $\sigma' = \sigma + \alpha p$（孔隙压力使有效拉应力增加）；
   压正 $\sigma' = \sigma - \alpha p$（教科书里最常见的形式）。

---

## 6. 指标记法与两个恒等式

采用爱因斯坦求和约定：重复指标自动求和。

$$
\delta_{ij} =
\begin{cases}
1, & i=j\\
0, & i\neq j
\end{cases}
$$

常用恒等式：

$$
\sigma_{ij}\delta_{ij} = \sigma_{ii} = \sigma_{11}+\sigma_{22}+\sigma_{33}
$$

$$
n_i n_i = 1,\qquad \frac{\partial n_i}{\partial n_j} = \delta_{ij}
$$

---

## 7. 数学工具清单

读本库只需要下面这些，遇到不熟的可以跳回来看：

- 矢量点积与投影： $\mathbf{a}\cdot\mathbf{b} = a_ib_i$
- 矩阵乘法： $(\mathbf{A}\mathbf{b})\_i = A_{ij}b_j$
- 特征值问题： $\mathbf{A}\mathbf{x} = \lambda\mathbf{x}$，等价于 $\det(\mathbf{A}-\lambda\mathbf{I})=0$
- 实对称矩阵：特征值全为实数，不同特征值的特征向量正交
- 三角函数： $\cos 2\theta = \cos^2\theta - \sin^2\theta$， $\sin 2\theta = 2\sin\theta\cos\theta$

---

下一篇 → [面力与应力张量](../01-stress/01-traction-and-stress-tensor.md)
