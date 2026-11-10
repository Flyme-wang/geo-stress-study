# M2 应变与变形

> **状态：提纲（待写）**
> 本页给出学习目标、推导路线与计划图件。写作时按 [`writing-style.md`](../00-notation/writing-style.md) 的六节结构展开。

---

## 学习目标

1. 由位移场定义小应变张量，并说明为什么需要对称化（转动与变形要分开）。
2. 求主应变与主方向，理解"应变也是张量"。
3. 区分**体积应变**与**偏应变**，与应力的球/偏分解对应。
4. 说明为什么应变不能任意给定（变形协调条件）。
5. 建立"应力–应变"配对使用的意识：莫尔圆有应力的，也有应变的。

---

## 计划章节

| 序 | 标题 | 关键推导 | 图件 |
|---|---|---|---|
| 02-1 | 位移梯度与应变张量 | $u_i$ 的梯度分解为对称+反对称；$\varepsilon_{ij}=\frac12(\partial_j u_i+\partial_i u_j)$ | 正方形微元的剪切/转动对比 |
| 02-2 | 主应变与应变不变量 | 与应力完全同构的特征值问题 | 应变莫尔圆 |
| 02-3 | 体积应变与剪胀 | $\varepsilon_v=\varepsilon_{kk}$；岩石剪胀（dilatancy） | 密砂/松砂的体积曲线 |
| 02-4 | 变形协调 | 二维协调方程 $\partial^2\varepsilon_x/\partial y^2+\dots=0$ | — |

---

## 关键结论（先记住，推导待补）

- 小应变张量对称，独立分量 6 个；反对称部分对应刚体转动 $\omega_{ij}=\frac12(\partial_j u_i-\partial_i u_j)$。
- $\varepsilon_v=\varepsilon_{11}+\varepsilon_{22}+\varepsilon_{33}=\Delta V/V$（小变形）。
- 各向同性线弹性下，$p=K\varepsilon_v$，$s_{ij}=2Ge_{ij}$——与应力分解一一对应。
- 岩石在三轴压缩中常先压缩后扩容（剪胀），这是压裂与巷道变形的重要机制。

---

## 前置

[M1 应力状态](../01-stress/README.md)、[符号与约定](../00-notation/symbols-and-conventions.md)

## 后续

[M3 本构与强度](../03-constitutive/README.md)
