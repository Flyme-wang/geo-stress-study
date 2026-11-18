# 案例拆解：《完全拉格朗日格式的求解流程（一）：公式推导与离散化》

> 来源：微信公众号「仿真札记」，2026-11-16
>
> 篇幅：一篇公式推导长文、6 张技术插图
>
> 对应知识图谱节点：`C6`
>
> 用途：几何非线性有限元的课前导读；CPS3 从形函数到 Newton 切线的完整链路
>
> 本库正文：[M6-10 完全拉格朗日格式](../docs/06-numerics/total-lagrangian-cps3.md)

---

## 0. 结论先行

这篇推文做了一件教材常拆散来讲的事：**用一个三节点三角形，把 TL 格式从形函数一直推到 Newton 切线**。它的链条是完整的——参考构形、位移梯度、变形梯度、Green–Lagrange 应变、B 矩阵、St. Venant–Kirchhoff 应力、虚功内力、材料与几何切线、迭代流程——中间没有跳步，也没有用"显然"带过关键项。

逐项核对后，原文的公式没有发现数学错误。$\mathbf{B}=\mathbf{B}\_L+\mathbf{B}\_{NL}$ 的第三行、$\mathbf{K}\_{\mathrm{geo}}=V\_0\mathbf{G}^{\mathsf T}\mathbf{M}\mathbf{G}$ 的分块结构都与独立推导一致。它的价值在于把"几何非线性"从一句口号变成可逐行核对的矩阵。

本库补强的是三件事：**数值可证伪**（每个公式都有 patch 检验）、**口径可复查**（工程剪应变与张量剪应变、平面应力的 $E\_{ZZ}$）、**实现可复现**（参考代码与校验脚本）。另外，$\mathbf{B}$ 第三行最容易出的"线性项重复计数"错误，本库在实现时真实踩过一次，记在 §3.1，供后来者对照。

适用范围：原文与本文都限于二维平面应力、St. Venant–Kirchhoff 材料、死载荷、无接触。三维、塑性、接触与压力随动载荷不在本文范围。

---

## 1. 原文结构与论证链

| 阶段 | 原文讲什么 | 本库位置 |
|---|---|---|
| 一、参考构形 | TL 与 UL 的区别；为什么在初始构形上积分 | [M6-10 §1.1–§1.2](../docs/06-numerics/total-lagrangian-cps3.md) |
| 二、几何与形函数 | CPS3 的节点、面积坐标、常数导数 | [M6-10 §1.3](../docs/06-numerics/total-lagrangian-cps3.md) |
| 三、位移梯度 | 插值求导，$\theta=\mathbf{G}\mathbf{d}$ | [M6-10 §2.1](../docs/06-numerics/total-lagrangian-cps3.md) |
| 四、变形与应变 | $\mathbf{F}=\mathbf{I}+\partial\mathbf{u}/\partial X$；Green–Lagrange 应变 | [M6-10 §2.2](../docs/06-numerics/total-lagrangian-cps3.md) |
| 五、B 矩阵 | 逐行求导；线性部分与非线性部分的分解 | [M6-10 §2.3](../docs/06-numerics/total-lagrangian-cps3.md) |
| 六、材料 | 平面应力 St. Venant–Kirchhoff；应力–应变矩阵 | [M6-10 §2.4](../docs/06-numerics/total-lagrangian-cps3.md) |
| 七、虚功与内力 | 参考构形上的虚功；$\mathbf{f}\_{\mathrm{int}}=V\_0\mathbf{B}^{\mathsf T}\mathbf{S}$ | [M6-10 §2.5](../docs/06-numerics/total-lagrangian-cps3.md) |
| 八、切线 | 材料项与几何项；几何刚度的紧凑形式 | [M6-10 §2.6](../docs/06-numerics/total-lagrangian-cps3.md) |
| 九、求解流程 | 残差、切线、更新；Newton 迭代 | [M6-10 §2.7](../docs/06-numerics/total-lagrangian-cps3.md) |

这张表也回答了"读原文还是读本库"：原文适合第一遍跟推导，图 2、图 3、图 5 把公式之间的关系画得很清楚；要落笔写代码、要查某个符号的口径，用本库的对应小节与校验脚本。

---

## 2. 原文讲对的五处

### 2.1 先选描述框架，再谈本构

原文没有一上来写应力–应变，而是先说明 TL 的积分域是初始构形、未知量是当前位移、应力度量是第二类 Piola–Kirchhoff 应力。这个顺序是对的：几何非线性首先是一个**记账方式**的选择，材料模型是第二层。很多教材把两者混在一起讲，学生容易以为"用了 TL 就必须用复杂的材料"。

### 2.2 位移梯度是整条链的枢纽

原文把 $\theta=[u\_{,X},u\_{,Y},v\_{,X},v\_{,Y}]^{\mathsf T}$ 单独拿出来，写成 $\theta=\mathbf{G}\mathbf{d}$，并强调 $\mathbf{G}$ 只由参考几何决定。这一步是整条链的枢纽：往前接形函数，往后接 $\mathbf{F}$、$\mathbf{E}$、$\mathbf{B}$ 与 $\mathbf{K}\_{\mathrm{geo}}$。抓住这个 4×6 矩阵，后面的公式都有共同的入口。

### 2.3 Green–Lagrange 应变的二次项没有被省略

原文把 $E\_{XX}$、$E\_{YY}$ 与 $\gamma\_{XY}$ 的二次项完整写出，没有写成"小应变 + 高阶小量"。这是对的：几何非线性的全部内容就在这些二次项里。若在这里省略，后面再补几何刚度就接不上——几何刚度正是这些二次项对位移的二阶导数。

### 2.4 B 矩阵的分解与第三行

原文把 $\mathbf{B}$ 拆成 $\mathbf{B}\_L+\mathbf{B}\_{NL}$，并明确 $\mathbf{B}\_L$ 就是 CST 的应变矩阵。第三行的两项 $u\_{,X}N\_{a,Y}+u\_{,Y}N\_{a,X}$ 与 $v\_{,X}N\_{a,Y}+v\_{,Y}N\_{a,X}$ 也写全了。这两项是 $\gamma\_{XY}$ 中两个交叉乘积的导数，最容易漏；原文没有漏。

### 2.5 几何刚度的紧凑形式

原文给出 $\mathbf{K}\_{\mathrm{geo}}=V\_0\mathbf{G}^{\mathsf T}\mathbf{M}\mathbf{G}$，并用分块矩阵表示 $\mathbf{M}$。这个形式比逐项展开更接近实现：代码里只需要组装两个 $2\times 2$ 的应力块，再乘 $\mathbf{G}$，不必为每个节点对写 36 个表达式。原文的图 5 把材料项与几何项的并联关系画得很清楚。

---

## 3. 本库补强与校核的六处

### 3.1 B 矩阵第三行的重复计数：一个真实的实现陷阱

本库在实现 `b_matrix` 时踩过一次：先把线性项 $(N\_{a,Y},N\_{a,X})$ 填进第三行，随后又用 $(1+u\_{,X})N\_{a,Y}+u\_{,Y}N\_{a,X}$ 覆盖式相加，结果线性部分被算了两遍。症状很有迷惑性：程序不报错，$\mathbf{B}(0)$ 检查也过，但有限差分切线误差达到 0.475，数值 Jacobian 不对称，Newton 要 24–30 次才勉强收敛。

正确写法是：线性项只填一次，非线性项只加 $u\_{,X}N\_{a,Y}+u\_{,Y}N\_{a,X}$ 与 $v\_{,X}N\_{a,Y}+v\_{,Y}N\_{a,X}$，不要再带那个 $1$。改正后，切线有限差分误差降到 $2\times 10^{-11}$，一致切线在 6 次内收敛。这个坑说明：**几何非线性的 B 矩阵不能靠"看起来像"来写，必须用 $\mathbf{B}(0)=\mathbf{B}\_L$ 与有限差分两条独立检查卡住。**

### 3.2 工程剪应变与张量剪应变的配对

原文的推导以 $\gamma\_{XY}=2E\_{XY}$ 为第三分量。这是常见选择，但必须在材料矩阵里成对写对：若应变向量是 $[E\_{XX},E\_{YY},\gamma\_{XY}]^{\mathsf T}$，则 $\mathbf{C}$ 的第三对角元是 $G=E/[2(1+\nu)]$；若改成张量剪应变 $[E\_{XX},E\_{YY},E\_{XY}]^{\mathsf T}$，同一位置要写 $2G$。两种写法都自洽，混用则剪切刚度差一倍，而且不会报错。

本库的检查办法是取纯剪状态，直接比较 $S\_{XY}$ 与 $G\gamma\_{XY}$。这比盯着矩阵看更可靠。

### 3.3 平面应力不等于厚度应变为零

原文的平面应力矩阵是对的，但没有展开 $E\_{ZZ}$。由三维 St. Venant–Kirchhoff 关系令 $S\_{ZZ}=0$，得到

$$
E_{ZZ}=-\frac{\nu}{1-\nu}(E_{XX}+E_{YY})
$$

取 $\nu=0.3$、$E\_{XX}=0.22$、$E\_{YY}=0$，则 $E\_{ZZ}=-0.0943$。也就是说，厚度方向收缩了近 10 %。CPS3 用常数厚度，不跟踪这个变化；面内刚度不受影响，但体积、厚度与压力载荷会受影响。把平面应力当成"厚度不变"，会在压裂、压实与软材料问题里低估体积变化。

### 3.4 几何刚度的分量核对

原文给出了 $\mathbf{K}\_{\mathrm{geo}}$ 的紧凑形式，本库补上了分量核对：取 $u\_a$ 与 $u\_b$ 两个自由度，由 $\mathbf{B}$ 的三行对 $u\_b$ 求导，乘应力求和，得到

$$
K_{\mathrm{geo}}^{u_au_b}=V_0\left[
S_{XX}N_{a,X}N_{b,X}
+S_{YY}N_{a,Y}N_{b,Y}
+S_{XY}(N_{a,X}N_{b,Y}+N_{a,Y}N_{b,X})
\right]
$$

它正好等于 $\mathbf{G}^{\mathsf T}\mathbf{M}\mathbf{G}$ 的对应项。这个核对的价值在于：**实现时不必直接组装 4×4 的 $\mathbf{M}$，可以按节点对写 2×2 的应力块；两种写法必须给出同一个矩阵。**

### 3.5 一致切线与仅材料切线的差距是数量级

原文说明了切线要包含几何项，但没有给出收敛数字。本库的 patch 检验给出：

- 一致切线：6 次迭代，残差 $1.41\times 10^{-1}\to 3.13\times 10^{-17}$；
- 仅材料切线：30 次迭代，残差停在 $5.23\times 10^{-9}$；
- 载荷从 0.2 加到 0.5，仅材料切线发散，末残差 2.10。

这张表把"一致切线"从术语变成可观察的现象：残差是否按平方下降，一眼就能判断切线有没有漏项。

### 3.6 St. Venant–Kirchhoff 的适用范围

原文把它作为几何非线性的标准材料，这是教学上的合理选择：能量形式简单、切线对称、代码短。但它不是通用橡胶模型。拉伸超过 10 % 量级时，应力会非物理地持续增大；某些变形模式下切线还会失去正定性。本库在 §1.5 与 §4 明确写了这个边界：**St. Venant–Kirchhoff 适合小应变、大转动，不适合大拉伸。**

---

## 4. 口径统一

| 项目 | 本文口径 | 说明 |
|---|---|---|
| 坐标系 | 参考构形 $(X,Y)$，当前构形 $(x,y)$ | TL 的积分与导数都用 $X$；当前坐标只用于 $\mathbf{F}$ |
| 位移梯度 | $\theta=[u\_{,X},u\_{,Y},v\_{,X},v\_{,Y}]^{\mathsf T}$ | 四个分量按此顺序；不同教材可能按 $[u\_{,X},v\_{,X},u\_{,Y},v\_{,Y}]$ 排，矩阵会转置 |
| 应变向量 | $[E\_{XX},E\_{YY},\gamma\_{XY}]^{\mathsf T}$ | 第三项是工程剪应变；若用张量剪应变，材料矩阵要改 |
| 应力向量 | $[S\_{XX},S\_{YY},S\_{XY}]^{\mathsf T}$ | 第二类 Piola–Kirchhoff 应力；不是 Cauchy 应力 |
| 符号 | 拉为正，压为负 | 与全库一致 |
| 厚度 | 常数 $h\_0$ | CPS3 不更新厚度；$E\_{ZZ}$ 只作诊断 |
| 载荷 | 死载荷 | 压力随法向转动的问题需要载荷刚度，本文不展开 |
| 数值单位 | 长度 m、模量 Pa | 脚本内部统一用 Pa；正文表格按 GPa 或 MPa 显示 |

---

## 5. 数值核对表

全部数字由 `py -3 tools/verify_tl_cps3.py` 复算，末尾打印 `ALL CHECKS PASSED`。参考单元为 $(0,0)$、$(2,0)$、$(0,1)$，$A\_0=1$、$V\_0=0.5$。

| 检查项 | 结果 |
|---|---|
| $\mathbf{B}(0)$ 与 $\mathbf{B}\_L$ 的最大差 | 0 |
| 有限差分应变导数与解析 $\mathbf{B}$ 的最大差 | $1.19\times 10^{-10}$ |
| 30° 刚体转动的 GL 应变最大分量 | $2.78\times 10^{-17}$ |
| 30° 刚体转动的小应变 $E\_{XX}$ | $-0.133975$ |
| 均匀拉伸 $s=0.2$ 的 $E\_{XX}$ | 0.220000 |
| 均匀拉伸的 $S\_{XX}$ | 17.2820 GPa |
| 均匀拉伸的 $S\_{YY}/S\_{XX}$ | 0.33 |
| 虚功比值 | 1.000000 |
| 应变能密度 | 1901.0212 MJ/m³ |
| 切线有限差分相对误差（三组） | $1.63\times 10^{-10}$、$1.99\times 10^{-12}$、$2.96\times 10^{-11}$ |
| 切线最大不对称量 | 小于 $10^{-12}$ |
| 零应力时 $\mathbf{K}\_{\mathrm{geo}}$ | 0 |
| 一致切线 Newton | 6 次，末残差 $3.13\times 10^{-17}$ |
| 仅材料切线 Newton | 30 次，末残差 $5.23\times 10^{-9}$ |
| 载荷 0.5、仅材料切线 | 发散，末残差 2.10 |

这些数字是**本库的可复现实验**，不是从原文抄来的。原文是符号推导，没有给出这组 patch 数字；两边可以独立对照的是公式结构，不是数值。

---

## 6. 6 张原图清点：登记 6、重绘 4

原文的 6 张插图按出现顺序提取，全部是技术示意图，没有装饰图。仓库是公开的，**原图一张都不复制**，只在表里登记远程地址，供读者回原文对照。本库图 47–图 50 是独立重绘，数据来自本库脚本，不是描图。

| 编号 | 主题（按原文上下文与图序登记） | 远程地址 |
|---|---|---|
| 01 | CPS3 参考构形：三节点三角形、节点坐标与厚度 | [原图](https://mmbiz.qpic.cn/mmbiz_png/8lRJyUQic2YBMtHWicbzKXsFycN1UZlEicCBCR2pOxaCEE7yrpv3ZE0orq5Yiaqadu0PxrZPRDcqa4AWvMDBM7MjtEB2M4lhAjMfV0rlygw35OI/640?wx_fmt=png&from=appmsg&watermark=1) |
| 02 | 位移插值 → 位移梯度 → 变形梯度 → Green–Lagrange 应变 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_png/8lRJyUQic2YAXzoVznKL9iaPnguhfXIxLtHaIdLqsVn4G3YRKoqqNnluu4L2QODV0ZD4dic4D39rbwvhwJneuwSQRUQtN00mpS3ahF5QWChX9Q/640?wx_fmt=png&from=appmsg&watermark=1) |
| 03 | 应变–位移矩阵 $\mathbf{B}$ 的组装：$\mathbf{B}\_L+\mathbf{B}\_{NL}$ | [原图](https://mmbiz.qpic.cn/mmbiz_png/8lRJyUQic2YC2oGmj1RpFibkxH0qVrVZGnxNL8u6ShkddvibibuTeMvBOiau2CRjj1eROGlibYvqanNAGhesiaJxx4zSXibNpdnwzPhOVWVChqB1S98/640?wx_fmt=png&from=appmsg&watermark=1) |
| 04 | 虚功原理 → 单元内力 $\mathbf{f}\_{\mathrm{int}}=V\_0\mathbf{B}^{\mathsf T}\mathbf{S}$ | [原图](https://mmbiz.qpic.cn/mmbiz_png/8lRJyUQic2YAX7gfTcPSCWicKCDGlSibWfm3H1DZrCGR6GLrLibYWVdOHjAhfrWKIwvqiaYlE4yd2rxHicj7Qxc1Mr2j6kCVMasLx4MMTZXenIbAs/640?wx_fmt=png&from=appmsg&watermark=1) |
| 05 | 切线刚度分解：材料项 $\mathbf{K}\_{\mathrm{mat}}$ 与几何项 $\mathbf{K}\_{\mathrm{geo}}$ | [原图](https://mmbiz.qpic.cn/mmbiz_png/8lRJyUQic2YDIKZ4O0s5Jgbc5m4Kd2kPCysyRCiba3vbu84RgRrqbNG8cBibXeAV8FkXNk3N6Wc60xWLEA3lbdiaiapZKr9pSKCmfvHbpqTJCkDg/640?wx_fmt=png&from=appmsg&watermark=1) |
| 06 | Newton–Raphson 迭代流程：残差、切线、更新 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_png/8lRJyUQic2YBFm1zJccBny58XRTzUZueliae7a08F8dHpShWQbw2XSNJiaXbUSkSOQQlu9fzNp2ZYx7BwDyxvwVQIsnhn6XZ7V0nfv7lO6kBA8/640?wx_fmt=png&from=appmsg&watermark=1) |

四张重绘图与原图的对应关系：

| 本库图件 | 对应原图 | 画的是什么 |
|---|---|---|
| 图 47 | 01 | CPS3 参考构形、面积坐标与 $\mathbf{G}$ 矩阵 |
| 图 48 | 02 | 构形映射、GL 应变客观性、有限拉伸二次项 |
| 图 49 | 03 | $\mathbf{B}\_L$、$\mathbf{B}\_{NL}$ 与总 $\mathbf{B}$ 的热图 |
| 图 50 | 04、05、06 | 切线分解、有限差分核对、Newton 残差与 patch |

三条处理原则：

1. **原图不复用**——公开仓库涉及版权，且微信图片缩放后不可读；
2. **重绘用本库数值**——图 50 的残差与误差由 `tl_cps3_lib.py` 实算，色标与坐标可追溯；
3. **编号独立**——本库图号连续（47–50），与原文图号无关，正文按需引用。

---

## 7. 课堂用法

| 时长 | 讲什么 | 怎么讲 |
|---|---|---|
| 45 min | 原文的图 1–图 3 与 M6-10 §0–§2.3 | 从"6 个节点自由度怎么变成 4 个位移梯度"开场；只写 $\mathbf{G}$ 与 $\mathbf{B}$，不碰切线 |
| 90 min | M6-10 §2.4–§2.8 | 板书 St. Venant–Kirchhoff 的二次项与几何刚度；让学生手算一个节点对的 $\mathbf{K}\_{\mathrm{geo}}$ |
| 4 学时 | M6-10 全篇 + 校验脚本 | 每人跑一遍 `verify_tl_cps3.py`，再故意删掉 $\mathbf{K}\_{\mathrm{geo}}$，观察残差从平方下降变成线性 |

四个现场接口：

- **橡胶与软材料**：小应变、大转动，St. Venant–Kirchhoff 可以作为第一近似；大拉伸换 Neo-Hookean。
- **岩土卸荷与基坑回弹**：几何非线性常被忽略，但深部卸荷的转动可能不小；先用本文的客观性检验判断小应变是否够用。
- **结构屈曲**：$\mathbf{K}\_{\mathrm{geo}}$ 不是修正项，而是主导项；特征值屈曲就是解 $\det(\mathbf{K}\_{\mathrm{mat}}+\lambda\mathbf{K}\_{\mathrm{geo}})=0$。
- **压裂与孔压耦合**：本文的 TL 公式是力学骨架；孔压耦合还要加 Biot 项与 Soils 步，见 [M6 导读](../docs/06-numerics/README.md) 的计划 06-2。

课堂提醒：不要让学生只背 $\mathbf{B}=\mathbf{B}\_L+\mathbf{B}\_{NL}$。让他们用 $\mathbf{B}(0)=\mathbf{B}\_L$ 与有限差分各验一次，比抄十遍公式有用。

---

## 8. 与知识图谱的对应

| 节点 | 状态 | 本案例贡献了什么 | 还需补什么 |
|---|---|---|---|
| `M6` 数值方法 | ✅ | 几何非线性从形函数到 Newton 的完整链路 | —— |
| `M6-07` 变分近似 | ✅ | 参考构形上的虚功与能量核对 | —— |
| `M6-08` 单元类型 | ✅ | CPS3 的形函数与平面应力假设 | —— |
| `M6-10` 完全拉格朗日格式 | ✅ | 公式推导、参考实现、数值检验 | 三维 TL、有限应变塑性、压力随动载荷 |
| `C6` 本案例 | ✅ | 6 张原图登记与 4 张图重绘 | —— |

后续接口：三维 TL（C3D8/C3D10）、有限应变弹塑性（乘法分解）、孔压耦合 TL、接触与摩擦的 TL 形式。本文只讲单相、二维、超弹性的一步，目的是把矩阵写对；真正的工程问题还要接本构与接触。

---

## 9. 引用与版权

- 原文：《完全拉格朗日格式的求解流程（一）：公式推导与离散化》，微信公众号「仿真札记」，2026-11-16。
  链接：[原文](https://mp.weixin.qq.com/s/vI0594vlIQejX1JZdO1Ftw)
- 本文为原文的独立改写与补充：CPS3 的几何、形函数、$\mathbf{G}$ 与 $\mathbf{B}$ 的逐项推导、
  平面应力的 $E\_{ZZ}$、$\mathbf{K}\_{\mathrm{geo}}$ 的分量核对、Newton 数值实验与误差表均为本库新增。
- **未复制原文任何图片**。6 张原图仅在 §6 登记远程地址；本库 4 张图件由
  [`tools/make_figures.py`](../tools/make_figures.py) 生成，
  数值由 [`tools/verify_tl_cps3.py`](../tools/verify_tl_cps3.py) 复算。
- 原文版权归原作者；链接仅作引用与出处标注。