# M1 综合练习（10 题，含答案）

> 使用建议：先不看答案，用纸笔完整写出步骤。每题都标注了考查的章节，
> 做不出来时回到对应页面，而不是直接看答案。

---

## A 组｜面力与张量

**A1.**（§01）某点应力张量（拉正，MPa）

$$
\boldsymbol\sigma=\begin{bmatrix}
40 & 15 & 0\\
15 & -10 & 0\\
0 & 0 & 25
\end{bmatrix}
$$

求法向 $\mathbf{n}=\frac{1}{\sqrt2}(1,1,0)$ 的截面上的面力矢量 $\mathbf{t}$，以及 $|\mathbf{t}|$。

**A2.**（§01）承接 A1，求该截面上的正应力 $\sigma_n$ 与剪应力 $\tau$。

**A3.**（§01）已知 $\sigma_{11}=30,\sigma_{22}=-10,\sigma_{12}=20$，其余为 0。
写出关于 $x_3$ 轴的力矩平衡式，检验对称性条件是否已经满足。

---

## B 组｜变换与主应力

**B1.**（§02、§03）求应力张量

$$
\begin{bmatrix}
80 & 40\\
40 & 20
\end{bmatrix}\ \mathrm{MPa}
$$

的主应力、主方向角，并画出莫尔圆的关键点（圆心、半径、两个主应力）。

**B2.**（§03）三维状态

$$
\begin{bmatrix}
50 & 30 & 0\\
30 & -10 & 0\\
0 & 0 & 20
\end{bmatrix}\ \mathrm{MPa}
$$

求三个主应力并排序。

**B3.**（§02）单轴拉伸 100 MPa。求正应力等于剪应力（$\sigma_n=\tau>0$）的截面角度。

---

## C 组｜莫尔圆与等效应力

**C1.**（§04）三维主应力 $(100,40,-20)$ MPa。求：
(a) 三个莫尔圆的圆心与半径；(b) $\tau_{\max}$ 及其所在截面的正应力。

**C2.**（§05）对 C1 的状态求 $J_2$、$\sigma_\mathrm{v}$，并用 Tresca 准则的等效应力量（$2\tau_{\max}$）作比较。

**C3.**（§05）判断下列说法对错，并说明理由：
(a) 静水压力下 $J_2=0$，所以应力张量为零；
(b) 纯剪 50 MPa 与单轴 86.6 MPa 的 von Mises 等效应力相同；
(c) 等效应力最大的位置一定是裂缝起裂位置。

---

## D 组｜地质应用

**D1.**（§06）埋深 2500 m，$\bar\rho=2450\ \mathrm{kg/m^3}$，$p_p=24.5$ MPa（正常压力），$\alpha=0.85$。
已知 $S_{h\min}=45$ MPa，$S_{H\max}=75$ MPa。求：
(a) $\sigma_v$；(b) 有效主应力；(c) 判断断层类型；(d) 检验水平应力比是否在 $\mu=0.6$ 的摩擦上限内。

**D2.**（§06）解释：为什么在逆冲断层环境中，水力压裂裂缝倾向于**水平**？

---

## 参考答案

<details>
<summary>A 组</summary>

**A1.** $\mathbf{n}=(1,1,0)/\sqrt2$。

$$
t_1=\sigma_{11}n_1+\sigma_{12}n_2=\frac{40+15}{\sqrt2}=38.89
$$

$$
t_2=\sigma_{21}n_1+\sigma_{22}n_2=\frac{15-10}{\sqrt2}=3.54
$$

$$
t_3=0
$$

$$
|\mathbf{t}|=\sqrt{38.89^2+3.54^2}=39.05\ \mathrm{MPa}
$$

**A2.** $\sigma_n=n_1t_1+n_2t_2=\frac{38.89+3.54}{\sqrt2}=30.0$ MPa。

$$
\tau=\sqrt{|\mathbf{t}|^2-\sigma_n^2}=\sqrt{1525-900}=\sqrt{625}=25.0\ \mathrm{MPa}
$$

**A3.** 关于 $x_3$ 轴：$\sigma_{12}-\sigma_{21}=20-20=0$，对称条件已满足（题给 $\sigma_{12}=\sigma_{21}$）。
若给出 $\sigma_{12}\ne\sigma_{21}$，则力矩不平衡，说明数据有误（或存在体力偶）。
</details>

<details>
<summary>B 组</summary>

**B1.** $C=(80+20)/2=50$，$R=\sqrt{30^2+40^2}=50$。

$$
\sigma_1=100,\quad \sigma_2=0
$$

$$
\tan2\theta_p=\frac{2\times40}{80-20}=\frac{4}{3}\Rightarrow 2\theta_p=53.13^\circ,\ \theta_p=26.57^\circ
$$

莫尔圆：圆心 $(50,0)$，半径 50，与横轴交于 0 与 100（过原点，因为 $\sigma_2=0$）。

**B2.** 先解 $x_1$–$x_2$ 子块：$C=20$，$R=\sqrt{30^2+30^2}=42.43$ → $62.43,-22.43$。
第三方向已对角，$\lambda_3=20$。排序：

$$
\sigma_1=62.43,\quad \sigma_2=20,\quad \sigma_3=-22.43\ \mathrm{MPa}
$$

**B3.** 设 $\sigma_n=\tau$：

$$
50+50\cos2\theta=50\sin2\theta
$$

$$
\sin2\theta-\cos2\theta=1
\Rightarrow \sqrt2\sin(2\theta-45^\circ)=1
\Rightarrow 2\theta=90^\circ\ \text{或}\ 180^\circ
$$

$2\theta=180^\circ$（$\theta=90^\circ$）给出 $\sigma_n=\tau=0$，非"$\tau>0$"。
故 $\theta=45^\circ$：$\sigma_n=50,\tau=50$。
</details>

<details>
<summary>C 组</summary>

**C1.** 排序 $100\ge40\ge-20$。

| 圆 | 圆心 | 半径 |
|---|---|---|
| $\sigma_1$–$\sigma_2$ | 70 | 30 |
| $\sigma_2$–$\sigma_3$ | 10 | 30 |
| $\sigma_1$–$\sigma_3$ | 40 | 60 |

$\tau_{\max}=60$ MPa，出现在与 $\sigma_1$、$\sigma_3$ 成 $45^\circ$ 的截面上，
该截面 $\sigma_n=40$ MPa。

**C2.** $J_2=\frac16[(60)^2+(60)^2+(-120)^2]=3600\ \mathrm{MPa^2}$，
$\sigma_\mathrm v=\sqrt{3\times3600}=103.9$ MPa。
Tresca：$\sigma_{\mathrm{T}}=2\tau_{\max}=120$ MPa。
两者不同（比值 $120/103.9=1.155$），在 $\sigma_2\neq(\sigma_1+\sigma_3)/2$ 时必然出现差异。

**C3.**
(a) 错。$J_2=0$ 但球张量不为零；总应力为 $p\delta_{ij}$。
(b) 对。$\sigma_\mathrm v=\sqrt3\times50=86.6$ 与单轴 86.6 相同。
(c) 错。脆性起裂取决于最大主应力（拉）或应力强度因子；等效应力是畸变指标。
在缺口/裂纹尖端，两者热点常不重合。
</details>

<details>
<summary>D 组</summary>

**D1.**
(a) $\sigma_v=2450\times9.81\times2500=60.1$ MPa。
(b) $\alpha p_p=0.85\times24.5=20.83$ MPa。

$$
\sigma_v'=60.1-20.8=39.3,\quad
S_{H\max}'=75-20.8=54.2,\quad
S_{h\min}'=45-20.8=24.2
$$

(c) 排序：$S_{H\max}'=54.2>\sigma_v'=39.3>S_{h\min}'=24.2$，竖直方向为中间主应力 → **走滑断层**。
(d) 比值 $54.2/24.2=2.24<3.1$，在 $\mu=0.6$ 的摩擦上限内，结果自洽。

**D2.** 逆冲环境下 $\sigma_3=\sigma_v$（垂向是最小主应力）。
裂缝垂直于最小主应力张开 → 裂缝面水平。这也是逆冲区压裂形成"水平缝"、施工压力高的原因。
</details>

---

返回 → [M1 应力状态](README.md)
