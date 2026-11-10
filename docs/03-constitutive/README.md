# M3 本构与强度

> **状态：提纲（待写）**
> 这一模块回答："材料在给定应力下如何响应？什么时候破坏？"

---

## 学习目标

1. 写出广义胡克定律，说明 $E,\nu$ 与 $K,G$ 的换算，讲清 $K,G$ 的物理意义。
2. 区分平面应力与平面应变，能算出 $\sigma_z$ 或 $\varepsilon_z$。
3. 掌握四种强度准则的适用对象与几何形状：
   - Rankine（最大主应力）——脆性拉裂；
   - Tresca（最大剪应力）——延性金属；
   - von Mises（畸变能）——延性金属；
   - Mohr–Coulomb——岩石、土、混凝土；
   - Hoek–Brown——节理岩体（经验准则）。
4. 能在 $\tau$–$\sigma$ 平面上画破坏包络线，并由莫尔圆判断是否破坏。
5. 理解塑性流动法则、剪胀角、软化/硬化对压裂与开挖的影响。

---

## 计划章节

| 序 | 标题 | 关键内容 |
|---|---|---|
| 03-1 | 广义胡克定律 | 四阶刚度张量、各向同性简化、$K$/$G$ 分解 |
| 03-2 | 平面应力与平面应变 | $\sigma_z$ 与 $\varepsilon_z$ 的换算；何时用哪个 |
| 03-3 | 屈服准则 | Tresca / Mises 的 $\pi$ 平面对比；$J_2$-$J_3$ 与 Lode 角 |
| 03-4 | Mohr–Coulomb | $c,\varphi$ 的物理意义；三轴试验如何标定；与 Mises 的关系（$\varphi=0$ 退化） |
| 03-5 | Hoek–Brown 与岩体强度 | GSI、扰动因子 $D$、$\sigma_{ci}$ |
| 03-6 | 塑性流动与剪胀 | 关联/非关联流动、剪胀角 $\psi$、软化 |

---

## 关键结论（先记住）

- 各向同性弹性：$K=\dfrac{E}{3(1-2\nu)}$，$G=\dfrac{E}{2(1+\nu)}$。
- Mohr–Coulomb：$\tau=c+\sigma'\tan\varphi$；主应力形式
  $\sigma_1'=\sigma_3'\dfrac{1+\sin\varphi}{1-\sin\varphi}+2c\dfrac{\cos\varphi}{1-\sin\varphi}$。
- $\varphi=0$ 时 MC 退化为 Tresca（$\sigma_1-\sigma_3=2c$）。
- 岩石的 $c,\varphi$ 必须用**有效应力**标定；用总应力标定会把孔隙压力混进强度参数。

---

## 前置

[M1 应力状态](../01-stress/README.md)、[M2 应变](../02-strain/README.md)

## 后续

[M4 断裂与破坏](../04-fracture/README.md)、[M5 地质力学](../05-geomechanics/README.md)
