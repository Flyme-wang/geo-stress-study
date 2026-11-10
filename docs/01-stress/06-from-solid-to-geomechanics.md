# 从固体力学到地质力学：压正、地应力与有效应力

> 对应知识图谱节点：`M1-06`（连接 `M5` 的桥梁）
> 前置：[主应力](03-principal-stress.md)、[偏应力与等效应力](05-deviator-and-von-mises.md)
> 后续：[地应力场](../05-geomechanics/README.md)（提纲）

---

## 0. 结论先行

1. 岩石力学与地质学**通常取压为正**。同一组数，换约定后主应力数值反号、排序反转。
   90% 的"公式对不上"事故来自这里。
2. 地壳中的应力是**三向受压**：
   - 垂向应力 $\sigma_v$ 基本等于上覆岩层重量；
   - 两个水平主应力 $S_{H\max}\ge S_{h\min}$，由构造、埋深与岩石强度共同决定。
3. 起作用的不是总应力，而是**有效应力**：

$$
\sigma'_{ij}=\sigma_{ij}-\alpha\ p_p\ \delta_{ij}\quad(\text{压正})
$$

4. 断层的类型（正断、走滑、逆冲）由**哪个主应力是竖直的**决定（Anderson 分类）。
5. 岩石强度受摩擦控制，因此 $S_{H\max}/S_{h\min}$ 不能任意大——存在由摩擦角决定的**上限**（应力多边形）。

---

## 1. 前提与符号

| 符号 | 含义 | 说明 |
|---|---|---|
| $\sigma_v$ | 垂向总应力 | 压正，正值 |
| $S_{H\max},\ S_{h\min}$ | 最大/最小水平主应力 | 压正 |
| $p_p$ | 孔隙压力 | 压正 |
| $\alpha$ | Biot 系数 | $0\lt \alpha\le1$，砂岩常取 0.7–1.0 |
| $\mu$ | 摩擦系数 | $\mu=\tan\varphi$ |
| $z$ | 埋深 | m |

本节**全部采用压为正**，与前面几页相反。原因：这一节的公式来自岩石力学文献，
强行改成拉正会与教材公式处处差号。**换算一次，之后统一用压正。**

---

## 2. 推导

### 2.1 第一步：垂向应力

取深度 $z$ 处一个水平面元，上方岩柱重量除以面积：

$$
\sigma_v(z)=\frac{1}{A}\int_0^z \rho(\zeta)\ g\ A\ \mathrm d\zeta=\int_0^z\rho(\zeta)g\ \mathrm d\zeta
\approx \bar\rho\ g\ z
$$

**依据**：静力平衡（上覆岩柱自重由下伏岩层支撑）。

取 $\bar\rho=2500\ \mathrm{kg/m^3}$， $g=9.81$：

$$
\sigma_v\approx 0.0245\ z\ \ (\mathrm{MPa},\ z\ \text{以 m 计})
\approx 24.5\ \mathrm{MPa/km}
$$

> 实际沉积盆地常见 22–26 MPa/km；欠压实泥岩可达 28 MPa/km 以上（超压）。

### 2.2 第二步：有效应力

Terzaghi（饱和土）：

$$
\sigma'=\sigma-p_p
$$

Biot（多孔岩石，考虑骨架可压缩性）：

$$
\sigma'_{ij}=\sigma_{ij}-\alpha p_p\delta_{ij},
\qquad \alpha=1-\frac{K_d}{K_s}
$$

其中 $K_d$ 为排水体积模量， $K_s$ 为矿物颗粒体积模量。

**依据**：孔隙压力作用于孔隙壁，抵消掉一部分总应力； $\alpha$ 描述"有多少孔隙压力真正作用在骨架上"。
$\alpha=1$ 退化为 Terzaghi； $\alpha=0$ 表示孔隙压力完全不影响骨架（致密无孔）。

**有效应力决定**：变形、强度、渗透率、裂缝张开与否。
**总应力决定**：与外部载荷的平衡。

### 2.3 第三步：Anderson 断层分类

三个主应力中，**哪一个竖直**决定了断层类型：

| 断层类型 | 竖直主应力 | 水平主应力关系 | 典型构造 |
|---|---|---|---|
| 正断层（normal） | $\sigma_1=\sigma_v$ | $S_{H\max}\gt \sigma_v\gt S_{h\min}$ | 伸展盆地、裂谷 |
| 走滑断层（strike-slip） | $\sigma_2=\sigma_v$ | $S_{H\max}\gt \sigma_v\gt S_{h\min}$（两者都水平） | 转换断层 |
| 逆冲断层（thrust） | $\sigma_3=\sigma_v$ | $S_{H\max}\gt S_{h\min}\gt \sigma_v$ | 造山带前缘 |

**依据**：断层沿最大剪应力面滑动，而最大剪应力面与 $\sigma_1$、 $\sigma_3$ 成 $45^\circ$
（[04](04-mohr-circle.md) §2.5）。 $\sigma_1$ 竖直 → 断层面倾角约 $60^\circ$（正断层）；
$\sigma_3$ 竖直 → 断层面约 $30^\circ$（逆冲）。

> 图 6：三类断层与主应力方位、断层面倾角示意。
> 文件：`assets/figures/fig06-anderson-faults.png`

### 2.4 第四步：摩擦平衡给出的应力上限

对无黏聚力（ $c=0$）的摩擦材料，破坏时莫尔圆与包络线 $\tau=\sigma'\tan\varphi$ 相切：

$$
\frac{\sigma_1'-\sigma_3'}{2}=\frac{\sigma_1'+\sigma_3'}{2}\sin\varphi
$$

**依据**：圆的半径 = 圆心到包络线的垂直距离（相切条件）。

整理：

$$
\frac{\sigma_1'}{\sigma_3'}=\frac{1+\sin\varphi}{1-\sin\varphi}
=\left(\frac{1+\sin\varphi}{\cos\varphi}\right)^2
=\left(\sqrt{1+\mu^2}+\mu\right)^2
$$

**依据**：三角恒等式 $\dfrac{1+\sin\varphi}{\cos\varphi}=\sec\varphi+\tan\varphi$。

对 $\mu=0.6$（ $\varphi\approx31^\circ$）：上限比值 $\approx3.1$。
即 $S_{H\max}'/S_{h\min}'\lesssim3.1$。这个上限加上" $\sigma_v$ 已知"，围出的区域就是**应力多边形**。

**工程含义**：如果你反演出的水平应力比超过 3，而该处没有强胶结或逆冲构造，
那这个结果是**不可信**的——地壳早就剪断了。

### 2.5 第五步：与水力压裂的接口（一句话预告）

裂缝沿**垂直于最小主应力**的方向张开（因为张开方向做正功最小）。
所以：

- 裂缝面法向 ≈ $S_{h\min}$ 方向 → 由裂缝方位可反推 $S_{h\min}$ 方位；
- 闭合压力 ≈ $S_{h\min}$ → 由压裂曲线可估算 $S_{h\min}$ 大小。

这两条是把井下数据"翻译"成应力张量的主要途径，详见 `docs/04-fracture/`。

---

## 3. 数值检验

**算例**： $z=2000$ m， $\bar\rho=2400\ \mathrm{kg/m^3}$， $p_p=19.6$ MPa（正常压力梯度 9.8 kPa/m），
$\alpha=0.9$， $S_{h\min}=35$ MPa。

1. $\sigma_v=2400\times9.81\times2000=47.1$ MPa。
2. 有效垂向应力： $\sigma_v'=47.1-0.9\times19.6=47.1-17.6=29.5$ MPa。
3. 有效最小水平应力： $S_{h\min}'=35-17.6=17.4$ MPa。
4. 比值： $\sigma_v'/S_{h\min}'=1.70\lt 3.1$，在摩擦上限内，合理。
5. 因为 $\sigma_v=47.1\gt 35=S_{h\min}$，若 $S_{H\max}\lt 47.1$ 则为**正断层**环境。

---

## 4. 常见误解

**误解 1："压正和拉正只是符号习惯不同，不影响结果。"**

为什么错：主应力排序、莫尔圆左右、 $\sigma_1-\sigma_3$ 的具体值全部改变。
若在压正公式里代入拉正的数，会得到" $\sigma_1\lt \sigma_3$"这类荒谬结果。

**误解 2：" $\sigma_v$ 一定等于 $\sigma_1$。"**

为什么错：只有在正断层环境（伸展）下才成立。逆冲环境下 $\sigma_v=\sigma_3$（最小主应力）。

**误解 3："孔隙压力只影响流体，不影响岩石强度。"**

为什么错：有效应力原理表明孔隙压力直接降低有效围压，从而降低抗剪强度。
这是水库诱发地震、注水诱发断层滑动的基本机制。

**误解 4："有效应力就是总应力减去孔隙压力，与 $\alpha$ 无关。"**

为什么错： $\sigma'=\sigma-\alpha p$ 中的 $\alpha$ 对致密岩石显著小于 1
（例如 $K_d/K_s\approx0.5$ 时 $\alpha\approx0.5$），忽略会高估有效应力。

---

## 5. 自测

1. 深度 3000 m、 $\bar\rho=2450$，求 $\sigma_v$。
2. 同上， $p_p=29.4$ MPa， $\alpha=0.8$，求 $\sigma_v'$。
3. 某地 $S_{H\max}=90$， $S_{h\min}=40$， $\sigma_v=60$ MPa。判断断层类型。
4. 用 $\mu=0.6$ 的摩擦上限，判断上题的水平应力比是否合理。
5. 为什么水力压裂裂缝通常垂直于 $S_{h\min}$ 而不是 $S_{H\max}$？

<details>
<summary>参考答案</summary>

1. $\sigma_v=2450\times9.81\times3000=72.1$ MPa。
2. $\sigma_v'=72.1-0.8\times29.4=72.1-23.5=48.6$ MPa。
3. $\sigma_1=S_{H\max}=90$、 $\sigma_2=\sigma_v=60$、 $\sigma_3=S_{h\min}=40$。
   竖直方向是中间主应力 → **走滑断层**。
4. $S_{H\max}'/S_{h\min}'=(90-\alpha p_p)/(40-\alpha p_p)$。
   取 $\alpha p_p=20$： $70/20=3.5\gt 3.1$，略超上限，需要更大的孔隙压力或存在胶结强度才合理。
5. 裂缝张开需要克服法向压应力。法向为 $S_{h\min}$ 的面承受最小的压应力，
   张开做的功最小，因此裂缝沿该面扩展（垂直于 $S_{h\min}$）。
</details>

---

## 延伸阅读

- Anderson, E. M. (1951). *The Dynamics of Faulting*, 2nd ed.
- Zoback, M. D. (2007). *Reservoir Geomechanics*, Ch. 1–4.（应力多边形与应力状态）
- Jaeger, Cook & Zimmerman (2007). *Fundamentals of Rock Mechanics*, Ch. 10.

返回 → [应力模块索引](README.md)
