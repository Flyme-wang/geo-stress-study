# 案例拆解：《从不可压缩到可压缩：N–S 方程与求解流程有哪些改变？》

> 来源：微信公众号「自由度FreeDof」，2026-10-29
>
> 篇幅：约 2 万字正文、8 章加 5 节附录、98 个独立公式、30 张插图
>
> 对应知识图谱节点：`C3`
>
> 用途：石油系统流体力学专题的课前阅读；可压缩求解器的自检清单
>
> 本库正文：[M7 流体力学与 CFD](../docs/07-fluid/README.md)

---

## 0. 结论先行

这篇长文值得完整吸收，原因不在"公式多"，而在**它按程序的执行顺序组织物理**：

> 存什么 → 怎么恢复 → 边界给什么 → 面上交换什么 → 怎么更新 → 怎么迭代 → 误差怎么查

它把"不可压缩→可压缩"的差别落到了七个具体动作上，每一步都给出可手算的数。
2 万字里没有跳步，这在公众号技术文里很少见。

本库按 M7 的四篇重新编排：控制方程与封闭、状态与边界、面通量与时间推进、
Sod 激波管走读。原文的论证链完整保留，另补了三类内容：

1. **原文没有展开的适用条件**：二维三维的边界计数、平均值恢复的偏差阶、正性保持的前提；
2. **原文没有做的定量校核**：网格收敛阶、Rusanov 在接触处的耗散、扩散步长系数、低马赫步数倍数；
3. **石油系统的接口**：井筒气侵、管线水击、超临界 CO2 输送分别对应哪一篇的哪个公式。

原文的物理判断没有发现错误；本轮做的是口径统一与补算，清单见 §3、§4。

---

## 1. 原文结构与论证链

| 原文章节 | 讲什么 | 本库位置 |
|---|---|---|
| 0 回顾 | 从离散变量到流场求解的总流程 | [M7-01 §2.1](../docs/07-fluid/01-governing-equations.md) |
| 1.1–1.3 | 质量、动量、能量三个方程的压缩性修正 | [M7-01 §2.2–§2.5](../docs/07-fluid/01-governing-equations.md) |
| 1.4 | 封闭：5 个方程、7 个未知量、2 条热力学关系 | [M7-01 §2.6](../docs/07-fluid/01-governing-equations.md) |
| 1.5 | 激波管作为共同检验 | [M7-04 §1](../docs/07-fluid/04-sod-shocktube.md) |
| 2.1–2.2 | 原始量与守恒量；由守恒量恢复压力温度 | [M7-02 §2.1–§2.2](../docs/07-fluid/02-states-and-boundaries.md) |
| 3.1–3.4 | 初始场、特征线边界、总静参数、壁面 | [M7-02 §2.3–§2.8](../docs/07-fluid/02-states-and-boundaries.md) |
| 4.1–4.6 | 物理通量、面状态、法向通量、Rusanov、黏性通量 | [M7-03 §2.1–§2.5](../docs/07-fluid/03-fluxes-and-time-stepping.md) |
| 5.1–5.5 | 半离散更新、显式、恢复、隐式、正性与守恒 | [M7-03 §2.6–§2.10](../docs/07-fluid/03-fluxes-and-time-stepping.md) |
| 6.1–6.2 | 稳态迭代与双时间步 | [M7-03 §2.8](../docs/07-fluid/03-fluxes-and-time-stepping.md) |
| 7.1–7.3 | 密度基/压力基、低马赫刚性、其他空间框架 | [M7-03 §2.9](../docs/07-fluid/03-fluxes-and-time-stepping.md) |
| 8.1–8.7 | 一次完整走读：初始化→通量→更新→恢复→误差→守恒 | [M7-04 §2–§3](../docs/07-fluid/04-sod-shocktube.md) |
| 附录 1–5 | 熵关系、声学关系、局部外推、波形、激波跃变 | [M7-02 §2.6](../docs/07-fluid/02-states-and-boundaries.md)、[M7-04 §2.2](../docs/07-fluid/04-sod-shocktube.md) |

这张表也是"读原文还是读本库"的答案：原文适合通读建立全局感；
要落笔算数、要查某一步的依据，用本库的对应小节。

---

## 2. 原文讲对了、可以直接搬进课堂的五处

### 2.1 一次更新里，所有面通量必须用同一层状态

原文在显式更新一节强调：先把当前层的全部面通量算完，再统一汇总更新，
不能"算一个面、更新一个单元、再算下一个面"。这句话是显式格式正确性的分水岭。
边算边更新会把同一时间层的数据混进下一层，格式退化成某种 Gauss–Seidel 迭代，
守恒审计立刻就对不上。

### 2.2 边界给几个条件，由信息传播方向决定

原文用外法向把特征值分成"进"和"出"，再数列数。这个讲法比"入口给速度、出口给压力"
的口诀可靠得多：口诀在超声速出口就失效，特征线计数不会。
本库把它做成四类边界的对照表（图 26），并补了二维三维的重数。

### 2.3 全域动量的增量等于两端压差乘时间

这是全文最漂亮的一个核对。质量与能量没有穿过边界，动量却因为压力一直在交换：

$$
\mathcal P(t)-\mathcal P(0)=(p_L-p_R)t .
$$

学生第一次看到"边界没有质量交换，但有动量交换"时，往往以为程序错了。
这个式子把"守恒"从口号变成可手算的 0.18。

### 2.4 低马赫数慢在时间步，不在网格

原文明确区分了两件事：声学 CFL 卡的是时间步，不是网格尺度；
网格该由空间分辨率决定，不能为了放大步长去粗化网格。这个区分很重要，
现场工程师常把"算得慢"归因于网格太密，实际是声速把步长钉死了。

### 2.5 守恒变量存的是平均值，恢复值不是体积平均

原文指出：除以密度、减动能都是非线性操作，由单元平均恢复出的压力、温度
一般不等于真实体积平均。这句话很多教材不写，但它是高阶格式必须做重构的原因。
本库把它量化成二阶小量的偏差（§3.2）。

---

## 3. 本库补强的六处

### 3.1 二维、三维的边界计数要加剪切模态

原文的计数表限定在一维，并明确说明了这一点；课堂上若直接推广到三维就会少算。
二维多一重剪切波，三维多两重：亚声速入口在二维要 3 条、三维要 4 条外部条件，
而出口始终只要 1 条。剪切条件给来流方向，不是给压力。
完整表格在 [M7-02 §2.4](../docs/07-fluid/02-states-and-boundaries.md)。

### 3.2 平均值恢复的偏差可以估阶

真实压力含的是动量的均方值除以密度，计算时只能用平均动量的平方除以平均密度。
由 Cauchy–Schwarz 不等式，

$$
\overline{\frac{\lVert\mathbf m\rVert^2}{\rho}}
\ge \frac{\lVert\bar{\mathbf m}\rVert^2}{\bar\rho},
$$

所以直接恢复会**高估**压力，差额是亚网格动能，光滑区为 $O(h^2)$。
一阶格式可以接受；高阶格式或强压缩区必须用重构后的面状态。
另一个陷阱是恢复顺序：必须先减动能再乘比热比减一，否则量纲都不对。
取密度 1 kg/m³、速度 300 m/s、压力 100 kPa，动能密度占总能量 15%，
忘了减就会把压力高估 18%。

### 3.3 Rusanov 在接触间断处的耗散有多大

原文说 Rusanov "对接触间断有较明显的抹宽"。本库把它算成一个数：
隔膜处精确质量通量 0.395391，Rusanov 给 0.517657，高 31%。
这解释了为什么一阶格式的接触间断会越抹越宽，也说明为什么需要 HLLC 或 Roe。

### 3.4 网格收敛阶要实算

原文给出三组网格与精确解的对照图，但没有列误差表。本库用 200、400、800 个单元
算密度积分误差：0.016646、0.010890、0.006978，观测阶 0.61、0.64。
低于一阶不是程序错：间断处的误差主导，一阶格式在接触与激波附近的收敛阶
接近 0.5 而不是 1。这条要讲给学生：不要拿一个低于理论阶的数字去怀疑代码，
也不要拿它去宣称"格式是二阶的"。

### 3.5 显式扩散步长要带系数

原文给量级 $\Delta t\sim h^2/D$。一维三点显式扩散的稳定界是 $h^2/(2\nu)$，
热扩散同理把运动黏度换成热扩散率。系数随差分模板与维数变化，
但"加密一倍、步长降到四分之一"这个比例不变。
本库用它算了甲烷：黏性限制约 630 s，声学限制约 $8.6\times10^{-6}$ s，
差七个数量级——所以显式可压缩格式的步长几乎总由声速决定。

### 3.6 低马赫的代价要算成倍数

原文说低马赫时显式方法"可能需要许多小步"。本库算成：甲烷马赫数 0.0295，
对流穿越一个单元要 $6.67\times10^{-4}$ s，声学步长 $8.61\times10^{-6}$ s，
二者相差 77.5 倍；若模拟 1 s 物理时间，显式约需 11.6 万步，
而对流本身只需 1500 步。这就是压力基与预处理存在的理由，
也是"低马赫刚性"四个字的定量含义。

---

## 4. 三处口径统一

| 编号 | 原文口径 | 本库口径 | 为什么统一 |
|---|---|---|---|
| 1 | 热扩散参数取 $\kappa/(\rho c_v)$ | 温度方程改用 $\kappa/(\rho c_p)$ | 同一个物理量的两种写法：内能扩散配定容比热，温度扩散配定压比热；常物性下只差一个比热比 |
| 2 | 扩散步长量级 $h^2/D$（略去系数） | $h^2/(2\nu)$、 $h^2/(2\alpha)$ | 手算稳定界要对得上；量级结论不变 |
| 3 | 边界示例 $c=300$ m/s、 $u=100$ m/s | 密度 1、压力 1、速度 0.5 的无量纲组 | 原文示例是有量纲示意；本库需要一组满足 $p=\rho RT$（其中 $R=1$）的数，才能把四步构造逐项复算 |

第 1、2 条不是错误，是工程文献的常见省略；第 3 条是为了让每个数都能在脚本里对上。
三处都写进了对应笔记的数值检验，正文没有"照抄原文再改口"的痕迹。

---

## 5. 数值核对表（本库全部复算）

所有数字可用 `py -3 tools/verify_fluid_notes.py` 重算，运行末尾应打印 `ALL CHECKS PASSED`。

| 项目 | 关键数值 | 核对的是什么 |
|---|---|---|
| 两侧声速 | 1.183216、1.058301 | 状态方程与比热比的一致性 |
| 隔膜 Rusanov 通量 | 0.517657、0.550000、1.331118 | 共享通量的三个分量 |
| 精确质量通量 | 0.395391 | 数值扩散的高估量（+31%） |
| 两单元更新 | 密度 0.948234/0.176766；总能 2.366888/0.383112 | 通量差与守恒 |
| 更新后恢复 | 压力 0.946328/0.150954；温度 0.997990/0.853975 | 减动能的顺序 |
| 星区状态 | 压力 0.303130、速度 0.927453 | 压力函数求根 |
| 波位置（t=0.2） | 0.263357、0.485945、0.685491、0.850431 | 稀疏波、接触、激波 |
| 网格收敛 | 误差 0.016646→0.010890→0.006978；阶 0.61/0.64 | 间断主导的误差 |
| 守恒审计 | 质量 0.5625、动量 0.18、能量 1.375；残差 $2.2\times10^{-16}$ | 边界收支 |
| 亚声速出口 | 压力 0.924385、速度 0.566080、面质量通量 0.564989 | 四步构造与通量 |
| 甲烷低马赫 | 声速 507.8 m/s、马赫数 0.0295、77.5 倍 | 声学与对流的时间尺度 |

这张表的作用和 FEM 那篇一样：**每一个公式都配一个能手算的量**。
学生算完能自己对，不必等老师批改。

---

## 6. 30 张原图清点：登记 30、重绘 10

原文 30 张插图按出现顺序提取，逐张核对主题，全部是技术示意图
（流程、几何、剖面、对照），没有装饰图混入。

仓库是公开的，**原文图片一张都不复制**；只在表里登记远程地址，供读者回原文对照。
本库的 10 张图（图 24–图 33）是独立重绘：数据来自本库脚本，不是描图。

| 编号 | 原图主题 | 远程地址 | 本库处理 |
|---|---|---|---|
| 01 | 可压缩流动的状态更新与热力学恢复 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/0G97QW1EFugx0WqVFhBOapYAMN45v4uFAGyMzTUcd8vM1oqUz4q4DGSv3Ho6lknMYtvDI1oPRmYoWkkQdAQHOm6DYiaD62jpVvewkibTwUDY4/640?wx_fmt=webp&from=appmsg) | 图 25 的「状态→通量→更新→恢复」闭环；未复制 |
| 02 | 固定控制体中的质量收支与密度变化 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFugA6rm4WTdEYMElcNMFMcuJTTYlib4NibDVICicoku3v9CyZIE5FMUmHdVibafDgE9HpRwU0m3cOUpfCic9icXvCnLRC5S2kmq7GIjia0/640?wx_fmt=webp&from=appmsg) | 图 24 左栏的控制体 |
| 03 | 守恒方程与模型封闭关系 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFujVWu7BAkho7XW1uX6rOliaPqeicW6nPicP5gzic1j0sqkicmibmjdQ9fU0pFTT41fZCj1YJZkJz96Xb9pV2fUadJmWLiawjlAtlH9ibso/640?wx_fmt=webp&from=appmsg) | 图 25 的封闭依赖关系 |
| 04 | Stokes 假设与两类黏性作用 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuhcPP7XWgE3CjCkvd81Ofs706qyxqelIFlaltAHLyLjV9z7Olq8Fwkujt8SoLA0R4MDQqLjUVk56QibtQLg1bfHufYCbuic9AzRA/640?wx_fmt=webp&from=appmsg) | 第 1 篇 §2.4 的公式与纯膨胀算例；未重绘 |
| 05 | 隔膜移去前后的激波管 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuia4aGibVVHdcNc7VxhyWNkEnCiaDBXn5UNiaNVpnL4srwVQiaAJpNIegzyYQ63psvUvnA7BZovNIPuMQia1ibGSmLMCybtxMpFCsN96I/640?wx_fmt=webp&from=appmsg) | 图 32 的左右初值 |
| 06 | 五个守恒分量与控制体收支 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/0G97QW1EFujJPa2hJCVaYXBvibmvhrkJle3HLp8CLF4tlgfxbGa9q3tfwpS20Cs27gWpaVIDiahPNTib5DoON6k0fc9nnRicco35o0CD1rbF7zY/640?wx_fmt=webp&from=appmsg) | 图 24 右栏的五个分量 |
| 07 | 由守恒状态恢复压力和温度 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuiaicZKXyybicBHvTiamyFw6UsSybvksL5oRW8TaFlKOvcEwLHFSPux7CoDJOQHKm2sfjG0pQVia8ObqNB86Jo7eK7F5t1EMFLFk6fE/640?wx_fmt=webp&from=appmsg) | 第 2 篇 §1.3 的逐步验算 |
| 08 | 沿外法向判断边界信息的进出 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/0G97QW1EFuia5arAM8iakmVEw9IPagib7t2dZnJkSRgtEsHdXelTK7aYI1nzsrBicianA7lVGlBpxIRA9b1I8z4O4SCSuR43KuhnA3BrJy0QJMM0/640?wx_fmt=webp&from=appmsg) | 图 26 的边界特征线 |
| 09 | 静态参数与理想停滞状态 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/0G97QW1EFuiaaibGlOq0ngscTYK7ZXtZqsiasZ45ibjibteOayI6JCCy2V2CNqLbQbHyiaD7AmzZGO0346OiaicRk0wdA2QDc42SxRR6c6XpEQ6bqRA/640?wx_fmt=webp&from=appmsg) | 第 2 篇 §2.5 的总静关系表；未重绘 |
| 10 | 亚声速出口的外侧状态构造 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuiav9uVXJzRTQYYFbCR9qunW3pkVvs6JVkfkhJk35g2xXDeG6ibv5cPjDdTibko1fMsFCKUNRvpt7MYGUmpYamsdq6fuLWicdL51Fo/640?wx_fmt=webp&from=appmsg) | 图 27 的四步构造 |
| 11 | 边界状态怎样进入面通量 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFujMMPcm9aNaEwcfe9ShgUuBUG2wMBic5hiayXy3S3mTtp9ek2B9GNVnQ4tK3icZVuVoWZfXqUM69BgiagdSDiaX2nyY6DE972Kva6icM/640?wx_fmt=webp&from=appmsg) | 图 28 的面通量进面 |
| 12 | 固定壁面的运动条件与热条件 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/0G97QW1EFujPtwqfiaoMLLiblRbBJg94ULmrIPFMiaZ0Y8m0Y5al4bne54C6xUlMLIz6XM9fgPo26qz4ZWAxMicGNRLywRB5N37vVSXz9vNticMU/640?wx_fmt=webp&from=appmsg) | 第 2 篇 §2.7 的壁面通量；未重绘 |
| 13 | 能量输运与压力做功 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/0G97QW1EFujrLgKepQ5psBe2rX3DK5OOU0Z2v4XNSF57ia25niapOW00ksEqwrJy3a9tiaziaMI5V93Duf2LRvbPkca9ic1S3ENFmOicvyK9x8CG8/640?wx_fmt=webp&from=appmsg) | 第 1 篇 §2.5 的内能方程推导；未重绘 |
| 14 | 单元平均值与左右面状态 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFugfG9q7UwicqByt5WwXNlnXIsEqzU5ibMgPHibwNleeKrLSibMhFk15TKDibGqObB9U0Sr4qC8pSFIhlMrls2R9HfuTatQe8sKiabYDo/640?wx_fmt=webp&from=appmsg) | 图 28 的左右面状态 |
| 15 | 斜向速度怎样穿过一个面 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuiatSXbCFOwf8mYkHqkDJnfawzNialyCm3I4bm5IHc31c0urX3V1G84modnLcfOcGHOvHnCOS4mpI46hqVcKcfxYobkBhbegPTUM/640?wx_fmt=webp&from=appmsg) | 第 3 篇 §2.2 的法向通量公式；未重绘 |
| 16 | Rusanov 共享通量的计算 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuiaAeupXiayr7RTB3U1iaglCKUZmic1an1QhejbdpRhnJKWwVJ45fiafMp6tOL6QNwFhZrnf2LavavczfyBv9EEdEatek9u4WYTL6K0/640?wx_fmt=webp&from=appmsg) | 图 29 的 Rusanov 通量 |
| 17 | 从面梯度得到应力与热流 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuhnaBPvKJvXvINupEW1bjH6y58pqCgfTWXsiap3XibWlDdl6JRQPbo20mQq7ia4kmvYY2Komqt6toK8CZM3j66SC2kd9ibTodkic0qg/640?wx_fmt=webp&from=appmsg) | 第 3 篇 §2.5 的黏性通量；未重绘 |
| 18 | 一个共享通量连接两个控制体 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuh4uTkejoibicIldH0K9GAldfa1MGY3vUkujJceMEsOk7roCDibjWZuXybIcibu0S1icdQkRDulU3iaDAGviag1lrg5cEwKv53Z7yNxqQ/640?wx_fmt=webp&from=appmsg) | 图 28 的两单元更新 |
| 19 | 显式更新与隐式迭代求解 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuhQxlrN7gRltrkSc7vnO7HbFFeWQf8u8iaNmSa3by6b0ZwXAbvzpBs3nIm6voxcIc7U9TNvRNeDweTHZWZVM7YibouF5x7xBeIKE/640?wx_fmt=webp&from=appmsg) | 图 30 的显式与隐式对照 |
| 20 | 缩小重构偏离，保留平均值 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/0G97QW1EFuia58fhicNibhmiaF5icJAO8r4mnFKbpRfnc2EVokGCxibjGm2buc5KT9sBkgUC3FkMzynI8PblduiaGzRFaKqqw2DaEYic0MQon99vGu8/640?wx_fmt=webp&from=appmsg) | 第 3 篇 §2.10 的正性保持；未重绘 |
| 21 | 物理时间与内层伪时间 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFugzzQHXfO8CrM7t6QF8SFDDaxz43wiaY55epILcJlkcsibwLCcKQsOpiaciclqthrHGarHNaDrF1xSwGZVwjc5bHmJlusjOO8ReovM/640?wx_fmt=webp&from=appmsg) | 图 30 右栏的伪时间迭代 |
| 22 | 密度基与压力基怎样确定压力 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFujMXc2icmAoGl7SCgibGz31EnTvIm1iaKmW16L2NeY1wzRygjwHibpFbKyWpeNL3K8LGHQIUSIgZHVNeAPloo5K4ON031VgUG8E0WE/640?wx_fmt=webp&from=appmsg) | 第 3 篇 §2.9 的密度基与压力基对照；未重绘 |
| 23 | 低马赫数下流动与声学的时间尺度 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFujofzuXVj3E368lVKpyOvdUsgTybYVMySHLMGeDx8O1D34IELUC73FwKqf0ibJfSnsOlGHrsSBVHmFiaEmOvWaxc39THmMsPTbQI/640?wx_fmt=webp&from=appmsg) | 图 31 的两个时间尺度 |
| 24 | 三个面的通量怎样更新两个控制体 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuhBW5h4IQDo5Q0uvBOkhrWoVDec1e0KYKQK4oHp71npu1icpQkkj3hneSYBzS4n8Xk2tf8NmKnfQKibiaqgkruX8Lbt2VmfKZxPX0/640?wx_fmt=webp&from=appmsg) | 图 28 的三面通量账 |
| 25 | Sod 激波管的数值解与精确解 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFujqP0RLL9hJibP9cTpq5zp7g2HYe0NicIib0LGhqseFaHpX3zZibxLvWXicn5HhiafoVNtuJwk9SOGia6pibov8ibNBMLuBhefcTib17KibFk/640?wx_fmt=webp&from=appmsg) | 图 32 的数值解与精确解 |
| 26 | 两端压差怎样改变全域动量 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/0G97QW1EFuj0ylSmzUia1B6SNmrBIRWpcrbklODWZFdE9Cv7INFT8iadpxXAC9buQxEPicbX0dJibukF8Kp5A25km3HkfEA0M29CfFOIQ7jnOMA/640?wx_fmt=webp&from=appmsg) | 第 4 篇 §3.6 的守恒审计表；未重绘 |
| 27 | 固定位置与跟随微团观察密度 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFuhENMkgA4aCplp7XWjMWia4ib9piauzSMJsia3p0CubBjXLXWOX1UtTqrsMcB0dBZxhnMF6ibQQgKDOtB7Mib6tjqlic4jBlH8ZxgLvgg/640?wx_fmt=webp&from=appmsg) | 第 1 篇 §2.2 的物质导数对比；未重绘 |
| 28 | 微团路径与声学信息路径 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/0G97QW1EFugiatus6akLWPc6lQPEIVn74uvasWJ8WRSibWFypLI6why6h6R3qo3ofiaIsadeh0cH63yY9qxrsAtZqolMNMg6lFhbsXibicibibUKm0/640?wx_fmt=webp&from=appmsg) | 第 2 篇 §2.6 的黎曼不变量；未重绘 |
| 29 | 沿声学路径追踪与同一时刻外推 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/0G97QW1EFugz9cdAa4kGKrfjlBgWWwnv5A67b5NXPdfGDj7Ql1paofZ0NQkXryx4kiaSkYy0TGuDFibY2juubDLBQHAYEOJYK2E92sqZDRxyQ/640?wx_fmt=webp&from=appmsg) | 图 27 的边界外推 |
| 30 | 随激波移动的控制体与质量收支 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/0G97QW1EFugBw7Fic7BmKn1YXCvwRntQLiciahaY1rqhLbVMZTAyONVTaVZVIRptyN2t5iclJKFJfjseb0e64UCcvslhhpIEoBXAICTibiafksnxg/640?wx_fmt=webp&from=appmsg) | 图 33 的激波控制体 |

十张重绘图与原图的对应关系：

| 本库图件 | 对应原图 | 画的是什么 |
|---|---|---|
| 图 24 | 02、06 | 控制体里存五个守恒量，面上交换三组通量 |
| 图 25 | 01、03、07 | 从守恒量到封闭方程组的依赖关系 |
| 图 26 | 08 | 外法向、特征线方向与进入条数 |
| 图 27 | 10、29 | 亚声速出口的四步构造 |
| 图 28 | 11、14、18、24 | 一个共享通量同时更新两个控制体 |
| 图 29 | 16 | Rusanov 通量：中心平均加数值扩散 |
| 图 30 | 19、21 | 显式推进与隐式迭代；物理时间与伪时间 |
| 图 31 | 23 | 低马赫数下声学与对流两个时间尺度 |
| 图 32 | 25 | Sod 激波管的数值解与精确解 |
| 图 33 | 30 | 随激波移动的控制体与质量收支 |

三条处理原则：

1. **原图不复用**——公开仓库涉及版权，且原文图为 WebP 截图，缩放后不可读；
2. **重绘用本库数值**——图 32 的三条曲线由 `tools/fluid_lib.py` 实算，色标与坐标可追溯；
3. **编号独立**——本库图号连续（24–33），与原文图号无关，正文按需引用。

---

## 7. 石油系统课堂用法

| 时长 | 讲什么 | 怎么讲 |
|---|---|---|
| 45 min | M7-01 与 M7-04 §0 | 从"压缩性什么时候重要"讲起，用马赫数 0.3 时 4.6% 的密度变化开场 |
| 90 min | M7-02 与 M7-03 §2 | 板书四类边界；让学生手算隔膜通量与两单元更新 |
| 4 学时 | M7 全模块 | 每人跑一遍校验脚本，再改 CFL 与通量格式看误差变化 |

三个现场接口：

- **井筒气侵与段塞**：多组分 Euler，验证顺序与 Sod 完全相同——先对精确解或解析解，
  再做守恒审计，最后看网格收敛；
- **管线水击**：声速决定波速与时间步，特征线边界就是 M7-02 的四类边界表；
- **超临界 CO2 输送**：状态方程换成非理想气体，恢复顺序与通量结构不变，
  但声速、比热必须由同一套热力学模型一致计算。

课堂提醒：不要把"可压缩"等同于"高速"。马赫数 0.3 时密度变化 4.6%，
已经不能忽略；但声速仍然决定时间步，这是两件不同的事。

---

## 8. 与知识图谱的对应

| 节点 | 状态 | 本案例贡献了什么 | 还需补什么 |
|---|---|---|---|
| `M7` 导读 | ✅ | 模块边界与石油接口 | —— |
| `M7-01` 控制方程与封闭 | ✅ | 5 方程 7 未知量的封闭账、不可压缩极限 | 非理想气体与真实气体状态方程 |
| `M7-02` 状态与边界 | ✅ | 恢复顺序、四类边界、四步出口构造 | 无反射边界与海绵层 |
| `M7-03` 面通量与时间推进 | ✅ | Rusanov、CFL、双时间、低马赫 | 重构与限制器、HLLC/Roe 的推导 |
| `M7-04` Sod 激波管 | ✅ | 精确解、误差阶、守恒审计 | 二维算例（双马赫反射等） |
| `C3` 本案例 | ✅ | 30 张原图的登记与 10 张图的重绘 | —— |

后续模块接口：M7-05 湍流与 RANS/LES、M7-06 多相与多组分已在导读中列出提纲；
储层渗流（达西）不在此模块，见 M5 与 M6-02 的孔压耦合。

---

## 9. 引用与版权

- 原文：《从不可压缩到可压缩：N–S 方程与求解流程有哪些改变？》，
  微信公众号「自由度FreeDof」，2026-10-29。
  链接：[原文](https://mp.weixin.qq.com/s/NZzyYyHdDG00ElneZVMAUA)
- 本文为原文的独立改写与补充：推导重排、数值算例、误差表、图件与练习题均为本库新增；
  原文措辞仅在转述观点时保留原意，未逐句复制。
- **未复制原文任何图片**。30 张原图仅在 §6 登记远程地址；本库 10 张图件由
  [`tools/make_figures.py`](../tools/make_figures.py) 生成，
  数值由 [`tools/verify_fluid_notes.py`](../tools/verify_fluid_notes.py) 复算。
- 原文版权归原作者；链接仅作引用与出处标注。
