# 案例拆解：《有限元单元类型》

> 来源：微信公众号「嘎子科研助手」，2026-11-14
>
> 篇幅：约 8 780 字、6 章、35 张插图
>
> 对应知识图谱节点：`C4`
>
> 用途：石油系统有限元课前阅读；单元选型与网格自检清单
>
> 本库正文：[M6-08 单元类型](../docs/06-numerics/element-types.md)

---

## 0. 结论先行

这篇推文没有公式，但它把一件很多教材讲不清的事讲清了：**单元类型是几何与变形模式的匹配问题**。作者从一辆自行车出发，先按"实心、薄、长、细"把零件分成实体、壳、梁、杆，再用一根悬臂梁解释"为什么厚度方向至少切四层"，最后用平面应力、平面应变、轴对称说明"少建一维"的三种约束。

它的论证链完整，工程判据可以直接搬进课堂：

- 薄：厚度小于面内特征长度的十分之一，用壳；
- 长：截面尺寸小于长度的十分之一，用梁；
- 细且只受拉：用杆，但受压必须另算稳定；
- 一阶实体在弯曲中：厚度方向至少 4 层；
- 复杂几何只能用四面体时：用二次四面体。

它没有展开的部分，正是本库补强的重点：剪切锁死的定量式、沙漏的秩与零能模式计数、壳/梁/杆的自由度与运动学假设、二维降维的第三维应力状态，以及原文那组 Abaqus 实测数字的口径。

原文的物理判断没有发现错误；需要提醒的是，**691 mm、2 层大三成、4 层大 7 % 是特定 Abaqus 设置下的实测**，不能当成普适常数。本库把可复现的部分与不可复现的部分分开写，见 §4、§5。

---

## 1. 原文结构与论证链

| 原文章节 | 讲什么 | 本库位置 |
|---|---|---|
| 引言 | 自行车与单元下拉菜单：C3D8R、C3D10、S4R、B31、T3D2 | [M6-08 §1.2](../docs/06-numerics/element-types.md) |
| 一、实体 | 曲柄切成六面体；节点三个平动；复杂形状用四面体 | [M6-08 §2.2](../docs/06-numerics/element-types.md) |
| 二、壳 | 挡泥板 1 mm × 600 mm × 100 mm；实体 24 万块 vs 壳 2400 块；中面网格加厚度输入 | [M6-08 §2.3](../docs/06-numerics/element-types.md) |
| 三、梁 | 车架钢管看成一维线；截面性质输入；焊缝细节丢失 | [M6-08 §2.4](../docs/06-numerics/element-types.md) |
| 四、杆 | 辐条 2 mm × 290 mm；约 2 kg 就弯，拉力可吊 100 kg；节点不转 | [M6-08 §2.5](../docs/06-numerics/element-types.md) |
| 五、同是方块分好几种 | 悬臂梁 100 × 10 mm；Q4 只有 2/3；剪切锁死；减缩积分；沙漏；4 层规则；Q8 99 %；640 块；三角形 23 %；C3D10 | [M6-08 §2.6–§2.7、§3.2–§3.5](../docs/06-numerics/element-types.md) |
| 六、少一维 | 链片平面应力；隧道平面应变；打气筒轴对称 | [M6-08 §2.8](../docs/06-numerics/element-types.md) |
| 结尾 | 自行车选型总结：曲柄实体、挡泥板壳、车架梁、辐条杆、链片平面应力 | [M6-08 §2.9](../docs/06-numerics/element-types.md) |

这张表也是"读原文还是读本库"的答案：原文适合第一遍建立直觉，尤其是挡泥板和辐条两个例子；要落笔算数、要查某个判据的适用范围，用本库的对应小节。

---

## 2. 原文讲对的五处

### 2.1 用"能不能表达变形"来选单元，而不是背名字

原文没有从定义出发，而是先问"这个零件主要怎么变形"：挡泥板主要受弯，所以厚度方向必须有足够的层；车架主要受弯和拉压，所以可以压成一根线；辐条只受拉，所以连转动都不需要。这个顺序比"先记单元名"可靠得多，因为名字会变，变形模式不会变。

### 2.2 挡泥板的自由度账

原文用 600 × 100 × 4 = 240 000 个实体块，对比 5 mm 中面网格的 2400 个壳单元，得出"少了一百倍"。这个账可以再往下算一层：实体 910 515 个自由度，壳 15 246 个，相差 59.7 倍。单元数差 100 倍、自由度差 59.7 倍，两个数字不矛盾，因为壳每个节点多两个转动自由度。

### 2.3 "直边弯不了，只能歪成梯形"

这是全文最漂亮的一句解释。它把剪切锁死从"数值现象"还原成几何事实：梁弯的时候，单元本该弯成弧，直边只能变梯形；梯形里量到的斜向错开，就是寄生剪切。学生一旦接受这个图像，后面"减缩积分""沙漏""二次单元"三个概念就能串成一条线。

### 2.4 减缩积分与沙漏的因果链

原文的链条是：直边 → 寄生剪切 → 只在中心点量 → 假变形量不到 → 锁死消失 → 但中心点也量不到真实弯曲 → 沙漏。这个顺序是对的，而且比"减缩积分精度高"这种口号有用。它解释了为什么减缩积分必须配沙漏控制，也解释了为什么"一层单元"在减缩积分下会错得离谱。

### 2.5 二次单元的"边能弯"

原文用"每条边中间再加一个节点，三点定一条弧"解释二次单元，这句话既准确又形象。它直接说明了 Q8 一层就能到 99 % 的原因：形函数空间里本来就有二次项，不需要靠多层网格去逼近弯曲。后面接"普通方块要 8 层 × 80 列 = 640 块"，对比也很干净。

---

## 3. 本库补强的六处

### 3.1 单元名字的解码

原文列出了 C3D8R、C3D10、S4R、B31、T3D2，但没有逐段解释。本库把它拆成"前缀 + 维度 + 节点数 + 后缀"：C 连续体、S 壳、B 梁、T 桁架；3D / 2D 是维度；数字是节点数；R 是减缩积分；I 是不协调模态；P 是孔压。这样遇到 C3D8P、C3D8I、S8R 这些名字，不需要逐个查表。

### 3.2 剪切锁死的定量式

原文说"直边只能歪成梯形"，本库把梯形里的剪应变算出来：纯弯曲下 Q4 的寄生剪应变是 $\gamma=\kappa a/(2\sqrt{3})$，与厚度无关。长厚比 10 时它是同点弯曲应变的 10 倍；按能量密度粗估，剪切项约为弯曲项的 38 倍。这个数解释了"为什么单元会显得硬"，也解释了"为什么长厚比越大越糟"。

### 3.3 沙漏的秩与零能模式计数

原文说"中心点什么变形都量不到"，本库把它写成秩的账：二维 Q4 有 8 个自由度，中心点只有 3 个独立应变分量，所以零空间维数是 5，其中 3 个刚体、2 个沙漏；三维 C3D8R 是 24 − 6 = 18，其中 6 个刚体、12 个沙漏。这个账可以直接在单刚特征值上验证（§5），也让"沙漏控制"不再是黑箱。

### 3.4 壳、梁、杆的自由度与运动学假设

原文说壳节点"除了挪还能转"，本库补全为：壳在整体坐标写 6 个自由度，公式用 5 个（3 平动加 2 个面内转动），绕法线的钻转由内部约束处理；梁每个节点 6 个自由度；杆只有 3 个平动。壳的厚度方向默认 5 个截面点，梁输入 A、I、J，杆只有 EA/L。这些细节决定了"什么时候壳够、什么时候必须实体"。

### 3.5 691 mm 的口径与 4 层规则的适用范围

原文的 691 mm、2 层大三成、4 层大 7 % 来自 Abaqus 默认设置下的减缩积分实体。本库不复现这三个绝对值，而是复现同一机制的可控量：

- 全积分 Q4 10×1 是 67.9 %，10×8 是 72.2 %；
- 单层 Q4 只沿长度加密到 160 块，也只到 91.5 %；
- Q8 10×1 是 99.8 %；
- T3 10×1 是 23.1 %；
- 640 块 Q4 是 99.8 %。

同时明确 4 层规则的适用范围：一阶、弯曲主导、减缩积分实体。全积分 Q4 10×4 也只有 72.0 %，说明"层数"不是万能的；二次单元一层就够。

### 3.6 二维降维的第三维应力状态

原文说平面应力"厚度方向没有应力"，平面应变"沿隧道方向不能伸缩"，这是对的，但没有给出第三维的应力量。本库补上：平面应变下法向应力等于泊松比乘以面内正应力之和；取两个面内正应力都是 −10 MPa、泊松比 0.3，法向应力是 −6 MPa；同样条件下平面应力的法向应力为零，法向应变是 $2.857\times10^{-5}$。这组数让学生看到"两种二维模型不是差一个系数，而是两个应力状态"。

---

## 4. 口径统一

| 编号 | 原文口径 | 本库口径 | 为什么统一 |
|---|---|---|---|
| 1 | 挡泥板 60 cm × 10 cm × 1 mm | 600 mm × 100 mm × 1 mm | 全文统一用 mm、N、MPa，避免单位混算 |
| 2 | "切一层""切四层" | 厚度方向的单元层数 | "层"只指厚度方向；面内加密不叫层 |
| 3 | "只有真值的三分之二""99 %""23 %" | 端部挠度除以 Euler–Bernoulli 真值 1.904762 mm | 原文没有写分母；本库固定分母，才能横向比较。二维平面应力真值比 Euler–Bernoulli 高约 0.8 %，不影响三位有效数字 |
| 4 | 691 mm、2 层大三成、4 层大 7 % | 标为"Abaqus 默认设置实测"，不逐位复现；另给全积分 Q4、Q8、T3 的可复现数 | 这三个数依赖网格、沙漏控制、加载方式与单元选项；直接当普适常数会误导 |
| 5 | "不到 2 公斤就弯""拉力吊 100 公斤" | Euler 屈曲 19.36 N，折合 1.97 kg；100 kg 拉力对应 312.2 MPa | 把手感换成可复算的临界力与应力 |
| 6 | "厚度不到宽度的十分之一" | 厚度小于面内特征长度的十分之一 | 挡泥板的"宽度"是 100 mm，但更一般的判据用面内特征长度 |

第 4 条不是错误，是实验口径。原文给的是"软件默认设置下会发生什么"，本库给的是"机制是什么、哪些量可以手算复现"。两者分工不同，课堂上要讲清楚。

---

## 5. 数值核对表

所有数字可用 `py -3 tools/verify_element_notes.py` 重算，运行末尾应打印 `ALL CHECKS PASSED`。

| 项目 | 关键数值 | 核对的是什么 |
|---|---|---|
| 梁真值 | 1.904762 mm | Euler–Bernoulli 基准 |
| Q4 10×1 / 10×2 / 10×4 / 10×8 | 0.678889 / 0.711090 / 0.720148 / 0.722471 | 只加层、不改长度 |
| Q4 80×8（640 块） | 0.997958 | 原文"8 层 × 80 列" |
| Q8 10×1 | 0.998259 | 原文"99 %" |
| T3 10×1 | 0.231178 | 原文"23 %" |
| 单层 Q4 20×1 / 40×1 / 160×1 | 0.842759 / 0.896881 / 0.915249 | 只沿长度加密的天花板 |
| 寄生剪应变 | 2.886751 | $\gamma=\kappa a/(2\sqrt{3})$ |
| 剪应变 / 弯曲应变 | 10.000000 | 长厚比 10 的锁死强度 |
| 零特征值个数 | 3 / 5 / 3 / 4 | Q4 全积分、Q4 减缩、Q8 全积分、Q8 减缩 |
| 最小正特征值 | 103 846 / 161 538 / 35 291 / 63 346 | 积分方案改变谱 |
| 实体 / 壳自由度 | 910 515 / 15 246 | 挡泥板，比值 59.72 |
| 辐条屈曲 / 应力 | 19.36 N（1.97 kg）/ 312.2 MPa | 2 mm × 290 mm 钢丝 |
| 平面应变 / 平面应力 | −6.0 MPa / $2.857\times10^{-5}$ | 第三维状态 |

这张表的作用和 FEM 那篇一样：**每一个结论都配一个能手算或能跑脚本核对的量**。学生不必等老师批改。

---

## 6. 35 张原图清点：登记 35、重绘 7

原文 35 张插图按出现顺序提取，逐张核对主题，全部是技术示意图、软件截图或场景图，没有装饰图混入。主题按原文上下文与图序登记，不转述图中文字。

仓库是公开的，**原文图片一张都不复制**；只在表里登记远程地址，供读者回原文对照。本库的 7 张图（图 34–图 40）是独立重绘：数据来自本库脚本，不是描图。

| 编号 | 主题（按原文上下文与图序登记） | 远程地址 |
|---|---|---|
| 01 | 单元类型下拉菜单：C3D8R、C3D10、S4R、B31、T3D2 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmQgr3uF9D3SLsq8SLjYAhCz3u9G8Njbd2tPee2eq8gicvoxBBAjibJ8jBwvyd4uWehb9VvxgYjcqGBDLfLOwILT8I2BHTyCfwN94/640?wx_fmt=jpeg&from=appmsg) |
| 02 | "全都切成小方块行不行"的场景图 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmS7ibkeeN1BGz3qPw6lbktMaj3VUY1aOLJ7H5ux34XCWFcib86GSYX7ZqSvQJIN6FSDicEzDtcmDSGppxiaUs3uryUen3ZsoZQAKlQ/640?wx_fmt=jpeg&from=appmsg) |
| 03 | 曲柄切成六面体，节点放在角上 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmS2NgxhU5vof7t5yqTiao6K0diaIvia2zpogsibsQFVYJ6uekia73icVibWORxpmsQUZlGtLb2JIPkqj8xwnyIEmeTxbXUh7hyvoXuKEs/640?wx_fmt=jpeg&from=appmsg) |
| 04 | 六面体节点的三个平动方向 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmSZcDTOEDLjrRygLJYpeZU9IwZ4WtD1lZ3W6lhdDdeoDw2PosGIwx82aSKRmicFyLPNByvMe6mGXclL4icgSQeN7ibYovL8WsMRqw/640?wx_fmt=jpeg&from=appmsg) |
| 05 | 实体单元：一块六面体 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmTicYmuFAPfnhpMVib3EZLc5wF5EdUn84nelWlfSKCmfSgxNWxAtctLZjWBedPiaicC5ibibjo4SicnKvNK7Au9uc5LCwc8vtIkty1IwM/640?wx_fmt=jpeg&from=appmsg) |
| 06 | 复杂形状改用四面体 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmReRy8RwPKnoItpziaFhAaLiatzBGU2mfhAyJ1rkMxvbrCTK4ibtONLZ8EkJ9w1ZKnBQt211423kHibeJ42kzHicianZ4VUoNrDVzErc/640?wx_fmt=jpeg&from=appmsg) |
| 07 | 挡泥板：1 mm 厚、600 mm 长、100 mm 宽 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmSWH6lmTW7Mo7diayUMdsM46nnBQSWAch7VjWOOPq1ibFM0eJfLHH8n2QlspZMEpCOsHPSmniayG1SibMoSKlBA6kukcgbNvFpRyVg/640?wx_fmt=jpeg&from=appmsg) |
| 08 | 厚度方向切四层：外侧拉长、内侧压短 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmTsblzEdoKFueBs0cgibIuzumfXibtBWJpdKUAzpfYe0nJiabeVib1KMUJsJL1jsJVHXuvwqjPQXrphaNkVXhIE6OWiaXkOFUCz9nFc/640?wx_fmt=jpeg&from=appmsg) |
| 09 | 600 × 100 × 4 = 24 万个方块 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmRbI62WuvOoiasNqsRPyBvUaMl2mgnG9v0FQPZFISEZibUOqEmF6dnSNHuKk75lEC9mmwDgZCsgJdIibiavicbZbsfOfJxGqnA7Dyqw/640?wx_fmt=jpeg&from=appmsg) |
| 10 | 中面铺 5 mm 网格：2400 个壳单元 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmTZiaLROtT0hwVw0vibHicA6WXbBu8x1cNyc5kpj2mwzX5jeJQT7ydWCgoicXHzMKHVkiae9NEib1u1W8sFMkZENWbG0aMXicBHfMicclc/640?wx_fmt=jpeg&from=appmsg) |
| 11 | 中面节点转动，厚度按输入值计算 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmTruu4K0DQooLGz45UCPiabSo7Gyd0ePTPAMVX4vFBIGpLX8libMWbxAJ0JySia7gaNUlic9InuXukWcYibL58J69rEbkxoq2yPYc6o/640?wx_fmt=jpeg&from=appmsg) |
| 12 | 壳单元 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmTZgcla6qokfBYzx8icMSTow7qvib5v4hkibzEeVib8cSmCAHQaqHFic5fOyer5OyI3S4sqvCiaIsrCD5EISBjB582s1xb5bVjEpRtxo/640?wx_fmt=jpeg&from=appmsg) |
| 13 | 车架钢管：600 mm 长、30–40 mm 粗，铺壳仍多 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmQWxDia7KGCcHjI5Yfkia9hhYsPPvqMkrutTRczFBDuLd22XfUEP66d8vvpvyicly8lfXW5WS5UczY3jByncOdrWokm1WPSXicfgjA/640?wx_fmt=jpeg&from=appmsg) |
| 14 | 梁：一根线加截面性质 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmTiccUwpnFf3wSrufibctkpicvEokIyoSiaVsTgGUBWIME7gQ5g64DH4XfwUN8V7l1Wmic5dAFMSRudFMTCPDs9s9U9b7wPia7esT9dg/640?wx_fmt=jpeg&from=appmsg) |
| 15 | 梁看不到焊缝细节，需局部改用实体或壳 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmRQceKicHiaOHMDOYqbeVuVicVtD2s3jhRIPSC6fGfsGUEGgTINO139DsNVgmaiawy05WUXu2B201AJJ1W6RyA7a2pCvicpaQSOsWuA/640?wx_fmt=jpeg&from=appmsg) |
| 16 | 辐条：2 mm 粗、290 mm 长 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmSKibftNibVTxmJAJvpGX3BybZbMGVYMJ8qHApCG1FZMXA5MJU81NvRWchsvEpLNRssEiago4OHfsfRazbhDiaxUp5ymwKoOJ14L0A/640?wx_fmt=jpeg&from=appmsg) |
| 17 | 杆单元：节点只平动，只算轴力 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmSXd0LEl6dMl9sicuMwHaKRWcPgqicuSUPexH7xXCltu9ic7j1hwaxvIFibSPgBUh4vedibjOxKsefVpwzs8RIE8rKiaAzc04jlHSKCw/640?wx_fmt=jpeg&from=appmsg) |
| 18 | 自行车上的分配：薄壳、长梁、细杆、实心实体 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmRAN9ia8z8EsA7icYyMaPecgegrrRACoUZP0axBA7MZLLq6XsgrBn3lTwN2KY7hAWnLaiahTGOAS3QIDiaLBvHicQxAbiawYNaFeVcL8/640?wx_fmt=jpeg&from=appmsg) |
| 19 | 悬臂梁基准：100 × 10 mm，一层 × 10 块 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmTLSbUkW5XOyYqzXIS7gPcbuvAgV6rJkcgP1UHgr6mof01nkibC1MO64CrzAJXLFuAkRUup9xKfKxeSd2EJpwpf6akw1ib6tW3LE/640?wx_fmt=jpeg&from=appmsg) |
| 20 | Q4 直边，端部只有真值的三分之二 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmQAr3dSn844a2cN4wWd8iam3jiadPllKG79SkdSS9aETRc7QiayEAGxpxhRyOYA0T05kZvrIeEqLwvwickNY4NAicJEq0xK8nLYcqcc/640?wx_fmt=jpeg&from=appmsg) |
| 21 | 直边弯成梯形：寄生剪切与剪切锁死 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmQnNssOtSfcX7IFeBCS8fCKDiaogFc1L74plP5sgIic1W6QTFKcyzzh6kkDdIzxEZz3K1MpgG13XyMWJpib1iaIHkvWXQjfN2W6bx0/640?wx_fmt=jpeg&from=appmsg) |
| 22 | 只在中心点积分：减缩积分 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmQRrmokVL20MXR3vWib5iaWoibibv52NoohNRF8EEiaUhanU1LMJdmDxmsXXoiaHKcNYJuLgO5EpN6DaSvpibe2yazNFrmLpas5YQUfgQ/640?wx_fmt=jpeg&from=appmsg) |
| 23 | 减缩积分后中线量不到变形：沙漏的起点 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmScJUAia8NQHiaHbe2y61UbG8hCb9I95jyPLxyE26M7Hu4VdEVpWrRkM3j8dErdvQOBYDtrGSIqP0KSm1CRaNjpGJhCG3EKsanEQ/640?wx_fmt=jpeg&from=appmsg) |
| 24 | 沙漏：网格交替上宽、下宽 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmQCfd8ia5Hg1J2XgHmLcdRSJaWT0eIicB9dHFrOXqL6uayJmbpCXRME2gM7ByK50ZMz4F03p4yocvRUryicav7IYA19ctd32P75K0/640?wx_fmt=jpeg&from=appmsg) |
| 25 | Abaqus 默认设置：1 层 691 mm、2 层大三成、4 层大 7 % | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmSqTic7AOTPezwn4nibWSvjHO2kicbrPlyXzegJlEHichq6HHMDwbYw7xIgymxUWy0aiczulOhsEhxQiapsvsiaRoqQ5NHFyJXjhK5IfA/640?wx_fmt=jpeg&from=appmsg) |
| 26 | 层数增加后每层中心离开中线，弯曲应变可测 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmRHlmvXskYjj2pm5ZXmpHrAwNrOoSuOTmz22rhtFZ9UbAhzE8J1p0Tq5TyAyLRibHaM6TRsvWU0crlz03VHeTb5RtnfrM2mzib04/640?wx_fmt=jpeg&from=appmsg) |
| 27 | 边中点加节点：二次单元 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmQdWEmXXY0adcCN8ibicNWU2AuQh30lFQre2CXOltnIDNnouhDOKIicL3uWYFcoGvUwe0jwtYxqibMUF6icuEdcTQeXvWp02BX1KDGo/640?wx_fmt=jpeg&from=appmsg) |
| 28 | Q8 一层 99 %；Q4 需 8 层 × 80 列 = 640 块 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmQkvicFiaR7bWOYia6jxmMm5hzmDGYgFEibkSiapdMLoq2fe9oPLb0es899SPTI7MEdmZFu2ibGUYeicMQskr7R1bhqVAgzgQubPJbkAM/640?wx_fmt=jpeg&from=appmsg) |
| 29 | 三角形一层只有 23 %；C3D10 是四角加六边中点 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmRWLpMl3eAHEpJV7CcEPu7v7gEBAdjOO88mRjlrDEZHGN1X5ibVnSLV88cPVKHQ8l58ObFiahU2ID6gdaspeHvzlMPQZwTkLzeYk/640?wx_fmt=jpeg&from=appmsg) |
| 30 | 平面应力：链片 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmRRtykyoKcpOlK14icJFd8slprnQOQcp7gw5ZWyBhU3yMYfLWkolb3k24okFujtW3e7Ku1xHDdTTibCjO98ms3FyiaBcUCib58YzBk/640?wx_fmt=jpeg&from=appmsg) |
| 31 | 平面应变：隧道截面 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmSzFSLsnlxicDGL07Vs7rzLcmd8xhQSiaEm7secjyxxsO4PyDZu1SgV2mxRNiaxBRCVWY4Lnaw1rafbgAxzot3dxSaKdpTBLxB9Hc/640?wx_fmt=jpeg&from=appmsg) |
| 32 | 轴对称：打气筒 | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmTmwQBv00TKWORgKTN4cf9Pu12QNWWHELFIvCU1icRUv50vpE2xa6IWqkOxqpuIJg5Y8gv0QTOEJ3EMGw4SbNsib3XIByFealCicw/640?wx_fmt=jpeg&from=appmsg) |
| 33 | 三种平面网格的差别在"没画出的方向" | [原图](https://mmbiz.qpic.cn/mmbiz_jpg/pBqL2KXaWmS8yIiaoPBEHPOYzwO3uYZLMibMicFgvUhT58MDAxxZVYYnOhia6OicsyeSqeJmb8lI61GcsSrqZLTEOSXzFmKf1gYj1gG0M2cfPWOs/640?wx_fmt=jpeg&from=appmsg) |
| 34 | 自行车选型总结：曲柄实体、挡泥板壳、车架梁、辐条杆、链片平面应力 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmRCU8Tsw7F5wVHrApZ1ib29omChsj6zaJLiaic4z8xDXNRhIGQSJ3R69NOIMY7OV6xreF7J1OAkq1TJQbnB3WZwku85ckTGLVPTtg/640?wx_fmt=jpeg&from=appmsg) |
| 35 | 下拉菜单的单元都是这几样：加节点、换量法、温度、孔隙水 | [原图](https://mmbiz.qpic.cn/sz_mmbiz_jpg/pBqL2KXaWmQJYxibibRib2F5MHyMwN7FlTZpgibF7oElOTK5ibeCiaImMiaLt0GEeX3VLvka7LrLCGyJzNYEuftsX0TjmyBphQJ6M8OACELibOK4oYo/640?wx_fmt=jpeg&from=appmsg) |

七张重绘图与原图的对应关系：

| 本库图件 | 对应原图 | 画的是什么 |
|---|---|---|
| 图 34 | 01、05、12、17 | 单元家族与 Abaqus 命名编码 |
| 图 35 | 03、04、08、11、15、17 | 实体、壳、梁、杆的几何与未知量 |
| 图 36 | 19、20、27、28、29 | 同一根梁、同样 10 块的端部精度 |
| 图 37 | 21 | 直边变梯形、寄生剪切与长厚比 |
| 图 38 | 22、23、24、26 | 减缩积分、零能模式与沙漏 |
| 图 39 | 30、31、32、33 | 平面应力、平面应变、轴对称 |
| 图 40 | 34、35 | 从几何到单元名的选型流程 |

三条处理原则：

1. **原图不复用**——公开仓库涉及版权，且原文图为微信压缩图，缩放后不可读；
2. **重绘用本库数值**——图 36 的挠曲线与图 37 的解析线由 `tools/element_lib.py` 实算，色标与坐标可追溯；
3. **编号独立**——本库图号连续（34–40），与原文图号无关，正文按需引用。

---

## 7. 石油系统课堂用法

| 时长 | 讲什么 | 怎么讲 |
|---|---|---|
| 45 min | 原文 §1–§4 与 M6-08 §0 | 从自行车四个零件开场：实心、薄、长、细分别对应什么单元；只讲几何与自由度，不写公式 |
| 90 min | M6-08 §2.2–§2.8 | 板书 Q4 的寄生剪切与减缩积分；让学生手算 2/3、10、5 个零特征值 |
| 4 学时 | M6-08 全篇 + 校验脚本 | 每人跑一遍 `verify_element_notes.py`，再改层数、积分方案与单元类型，看端部精度怎么变 |

四个现场接口：

- **套管与油管**：薄壁、内压加轴向力，用 S4 / S4R；长细段用 B31。内压产生环向应力，壳的截面点要够；局部螺纹与接箍必须拿实体子模型重算。
- **井架与桁架**：只受轴力，用 T3D2；受压杆必须另做屈曲分析。线性杆单元不会告诉你"这根杆已经失稳"。
- **盐构造与大变形**：用 CEL，欧拉区域用 EC3D8R，避免网格畸变；拉格朗日区域与欧拉区域用耦合面交换力。
- **水力压裂与储层压实**：裂缝用 COH3D8 / COH2D4，孔压用 C3D8P / CPE4P，配 Soils 步；裂缝开度是 cohesive 的损伤变量算出来的，不是"网格变形"直接给的。

课堂提醒：不要把"单元高级"当成"结果准确"。Q8 比 Q4 贵，但在接触问题里未必更稳；C3D10 比 C3D8 适应复杂几何，但在显式分析里时间步更小。选单元是精度、成本、稳定性三者的折中。

---

## 8. 与知识图谱的对应

| 节点 | 状态 | 本案例贡献了什么 | 还需补什么 |
|---|---|---|---|
| `M6` 数值方法 | ✅ | 单元选型与降维的工程判据 | —— |
| `M6-06` 雅可比矩阵 | ✅ | 等参映射在单元刚度里的位置 | —— |
| `M6-07` 变分近似 | ✅ | 离散解偏硬与剪切锁死的两条解释 | —— |
| `M6-08` 单元类型 | ✅ | 形函数、寄生剪切、沙漏、二维降维的完整推导 | 接触单元、XFEM 富集单元、CEL 的欧拉公式 |
| `C4` 本案例 | ✅ | 35 张原图的登记与 7 张图的重绘 | —— |

后续模块接口：孔压耦合与 Soils 步（06-2）、cohesive 与 XFEM（06-3）、CEL（06-4）已在 M6 计划中列出；本文只讲单元几何与积分方案，不展开本构与接触算法。

---

## 9. 引用与版权

- 原文：《有限元单元类型》，微信公众号「嘎子科研助手」，作者 北工大-嘎子，2026-11-14。
  链接：[原文](https://mp.weixin.qq.com/s/jhRFROb_Cc7wdrebnnD5ig)
- 本文为原文的独立改写与补充：形函数、应变矩阵、寄生剪切、沙漏、二维降维的推导与数值算例均为本库新增；
  原文的工程判据与 Abaqus 实测数字按原意转述，并在 §4 注明口径。
- **未复制原文任何图片**。35 张原图仅在 §6 登记远程地址；本库 7 张图件由
  [`tools/make_figures.py`](../tools/make_figures.py) 生成，
  数值由 [`tools/verify_element_notes.py`](../tools/verify_element_notes.py) 复算。
- 原文版权归原作者；链接仅作引用与出处标注。
