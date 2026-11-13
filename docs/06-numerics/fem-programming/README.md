# M6-P 有限元编程基础

> **状态：已完成（5 篇）**
> 对应知识图谱节点：`M6-P`（子节点 `M6-P1` ～ `M6-P5`）
> 前置：[M1 应力状态](../../01-stress/README.md)、[M3 本构与强度](../../03-constitutive/README.md)
> 后续：[M6-02 孔压耦合与 Soils 步](../README.md)、[M6-03 Cohesive 与 XFEM](../README.md)

这一组笔记回答一个具体问题：**把 $K u = F$ 写成程序，到底要写哪几步，每一步为什么这样写？**

它不依赖某个软件，也不假设读者学过 Fortran。所有算例都可以用
[`tools/verify_fem_notes.py`](../../../tools/verify_fem_notes.py) 重算：

```bash
py -3 tools/verify_fem_notes.py
```

---

## 0. 五篇的路线

| 章 | 标题 | 一句话结论 | 图 |
|---|---|---|---|
| [1](01-language-and-first-program.md) | 从语言到第一个程序 | 有限元程序只有四件事：离散、单元刚度、组装、求解；语言只影响写法，不影响数学 | 图 7、图 8 |
| [2](02-node-ordering-and-jacobian.md) | 节点顺序、形函数与雅可比 | 节点顺序决定"哪个形函数配哪个物理角点"；detJ 是面积缩放，必须处处为正 | 图 9、图 10 |
| [3](03-assembly-and-sparse-storage.md) | 组装与稀疏存储 | 组装是按自由度映射的"累加"，不是填表；稀疏存储省的是零，不是数 | 图 11、图 12 |
| [4](04-constraints-and-loads.md) | 约束与荷载 | 三种约束写法在精确算术下等价；三种荷载三个积分公式 | 图 13、图 14 |
| [5](05-linear-solvers.md) | 线性方程组求解 | 直接法稳、迭代法省内存；收敛看真实残差，不看"步数好看" | 图 15、图 16 |

学习顺序就是表的顺序。第 1 章建立全局图景，第 2、3 章是"单元 → 整体"的两个关键转换，
第 4 章把边界条件与右端项补齐，第 5 章才谈"怎么解、为什么解得起"。

---

## 1. 与原始推文的对应

本组笔记整理自微信公众号「有限元编程基础」系列 11 篇（2026-11），原文链接见
[案例笔记](../../../cases/2026-11-12-fem-programming-series.md)。对应关系：

| 原文 | 标题 | 本库位置 |
|---|---|---|
| 01 | 为什么选择 Fortran 语言？ | [第 1 章 §1.1、§4](01-language-and-first-program.md) |
| 02 | 有限单元法的本质 | [第 1 章 §2.1](01-language-and-first-program.md) |
| 03 | 从 1 个小型有限元程序开始 | [第 1 章 §2.4、§3.3](01-language-and-first-program.md) |
| 04 | 单元节点排序的重要性 | [第 2 章 §2.1、§3.1](02-node-ordering-and-jacobian.md) |
| 05 | 雅可比矩阵 | [第 2 章 §2.2–§2.4、§3.2](02-node-ordering-and-jacobian.md) |
| 06 | 刚度矩阵组装 | [第 3 章 §2.1–§2.3、§3.1](03-assembly-and-sparse-storage.md) |
| 07 | 刚度矩阵存储技术 | [第 3 章 §2.4–§2.6、§3.2](03-assembly-and-sparse-storage.md) |
| 08 | 固定约束处理 | [第 4 章 §2.1–§2.2](04-constraints-and-loads.md) |
| 09 | 位移约束处理 | [第 4 章 §2.2–§2.4](04-constraints-and-loads.md) |
| 10 | 形成荷载向量 | [第 4 章 §2.6–§2.7、§3.3](04-constraints-and-loads.md) |
| 11 | 线性方程组求解 | [第 5 章](05-linear-solvers.md) |

原文的长处是把 11 个主题串成一条完整的实现链；本库做的是**补依据、补数值核对、纠三处笔误**。
逐条对照与纠错清单见案例笔记，正文只在需要处引用原文结论。

### 1.1 本库对原文的三处订正

1. **自由度编号**：原文第 03 篇正文说"节点 3 的 $y$ 方向自由度为 2、节点 4 的为 4"，
   但同一篇给出的编号公式是 $\mathrm{dof}(c,i)=2(i-1)+c$，且荷载向量填在 6、8 位。
   正确值是 **6 和 8**。本库第 1 章 §2.3 用公式重新推一遍。
2. **CSR 的 PTR 长度**：原文说 PTR 长度等于自由度个数，示例给 `PTR(1:5) = [1,2,4,6,8]`；
   但要按 `PTR(i+1)-1` 取第 $i$ 行末位置，第 5 行会越界。完整 CSR 需补末位哨兵
   `PTR(6)=10`，长度 $n+1$。本库第 3 章 §2.6 与图 12 按完整写法给出。
3. **示例矩阵的对称性**：原文 5×5 示例存了 $K_{21}$ 却没存 $K_{12}$。
   存储格式本身允许非对称，但线弹性刚度矩阵是对称的，实际只需存一半。
   本库把它明确标为"格式示意，不是物理刚度矩阵"。

---

## 2. 统一约定

- **节点顺序**：逆时针，从左下角开始，记为 1→2→3→4。所有公式与图都按这个顺序。
- **自由度编号**：节点 $i$ 的 $x$ 方向为 $2i-1$， $y$ 方向为 $2i$。一维数组下标从 1 开始讲，
  Python 实现时整体减 1。
- **矩阵存储**：讲解用 1 基下标（与 Fortran、MATLAB 一致）；Python 代码用 0 基。
- **单位**：长度 m、力 N、应力 Pa；弹性模量给 300 GPa 时写作 $3.0\times10^{11}$ Pa。
- **算例材料**： $E=300$ GPa、 $\nu=0.3$、平面应力、厚度 $t$ 按题注明。
  这三个数在五篇里复用，便于交叉核对。

---

## 3. 图件

全部由 [`tools/make_figures.py`](../../../tools/make_figures.py) 生成：

| 图 | 文件 | 用在哪 |
|---|---|---|
| 图 7 | `assets/figures/fig07-single-element.png` | 第 1 章：单 Q4 的模型与 $u_x$、$\lvert u\rvert$、$\sigma_y$ 云图 |
| 图 8 | `assets/figures/fig08-fem-workflow.png` | 第 1 章：四步流程 |
| 图 9 | `assets/figures/fig09-node-ordering-detj.png` | 第 2 章：三种节点顺序与 detJ 符号 |
| 图 10 | `assets/figures/fig10-isoparametric-jacobian.png` | 第 2 章：母单元到实际单元的映射 |
| 图 11 | `assets/figures/fig11-assembly-scatter.png` | 第 3 章：两个 Q4 的组装与稀疏结构 |
| 图 12 | `assets/figures/fig12-sparse-storage.png` | 第 3 章：Skyline / COO / CSR |
| 图 13 | `assets/figures/fig13-constraint-methods.png` | 第 4 章：三种约束写法 |
| 图 14 | `assets/figures/fig14-load-vector.png` | 第 4 章：三类荷载、点荷载一致性与组装后的 F |
| 图 15 | `assets/figures/fig15-solver-cost.png` | 第 5 章：三种求解器的代价量级 |
| 图 16 | `assets/figures/fig16-cg-pcg-convergence.png` | 第 5 章：CG 与 PCG 的真实残差曲线 |

重绘：

```bash
py -3 tools/make_figures.py
```

---

## 4. 学完应该能做到

1. 拿到一个二维 Q4 网格，手写出它的自由度编号表与单元连接表。
2. 给定节点坐标，判断 detJ 是否处处为正；解释"负 detJ 不能取绝对值了事"。
3. 用两个单元的例子手算组装后的非零结构，并说明为什么共享节点的项要累加。
4. 分别用消元法、改行列法、罚函数法施加同一个位移约束，并说出三者的误差与条件数代价。
5. 把面力、体力、集中力分别写成一致的节点荷载，并用合力核对。
6. 解释为什么 10 万自由度以上要换迭代法，以及为什么 CG 的残差曲线不是单调下降的。

---

## 5. 延伸阅读

- Zienkiewicz, Taylor & Zhu, *The Finite Element Method: Its Basis and Fundamentals*, 7th ed.
  —— 弱形式与等参元的标准参考。
- Bathe, *Finite Element Procedures*, 2nd ed. —— 组装、带宽、求解器的工程视角。
- Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed. —— CG、PCG 与预处理。
- 原始推文系列（11 篇），链接与逐条拆解见
  [案例笔记](../../../cases/2026-11-12-fem-programming-series.md)。
