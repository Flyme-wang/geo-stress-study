# 参考插图（仅作对照）

本目录保存**原文引用的原始插图**，只用于核对我们自绘图件的几何关系，不参与正文推导，也不替代 `assets/figures/` 中的自绘图。

| 文件 | 内容 | 来源文献 | 许可 |
|---|---|---|---|
| `ref01-stress-components-ma17071498.jpg` | 三维微元三个截面上的正应力与切应力分量 | [10.3390/ma17071498](https://doi.org/10.3390/ma17071498) | CC BY 4.0（MDPI） |
| `ref02-mohr-circles-3d-ma17143436.jpg` | 三向主应力对应的三个莫尔圆与可行域 | [10.3390/ma17143436](https://doi.org/10.3390/ma17143436) | CC BY 4.0（MDPI） |
| `ref03-yield-criteria-applmech3010018.jpg` | 平面应力下 Rankine / von Mises / Tresca 的对比 | [10.3390/applmech3010018](https://doi.org/10.3390/applmech3010018) | CC BY 4.0（MDPI） |

## 使用规则

1. **不进入正文插图位**：正文只使用 `assets/figures/` 下由 [`tools/make_figures.py`](../../tools/make_figures.py) 生成的图，保证坐标系、符号与文字说明可控。
2. **不裁改、不加标注**：如需讲解某处细节，另绘一张图，不用修图软件在原图上涂改。
3. **引用时标注**：在幻灯片或讲义中使用时写清文献 DOI 与 “CC BY 4.0”。
4. **对照清单**：
   - `ref01` ↔ 自绘 [`fig01-stress-cube.png`](../figures/fig01-stress-cube.png)：核对三个可见面上的正应力方向与切应力下标。
   - `ref02` ↔ 自绘 [`fig03-mohr-circles-3d.png`](../figures/fig03-mohr-circles-3d.png)：原文图把三个圆心到原点的距离标成 $(\sigma_i+\sigma_j)/2$，容易与“圆心坐标”混淆；自绘图直接标 $C_{12}$、 $C_{23}$、 $C_{13}$。
   - `ref03` ↔ 自绘 [`fig04-yield-criteria.png`](../figures/fig04-yield-criteria.png)：原文图为归一化坐标 $\sigma_{a,1}/\sigma_{a,\mathrm{eq}}$，自绘图为 $\sigma/\sigma_y$，两者的形状相同、刻度含义不同。

## 版权说明

三张图均来自 MDPI 开放获取期刊，按 CC BY 4.0 许可转载，版权归原作者所有。本仓库不对其主张任何权利；正文中的推导、算例与自绘图件为独立完成。
