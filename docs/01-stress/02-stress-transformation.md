# 斜截面上的正应力与剪应力：100 MPa 为什么变成 50/50

> 对应知识图谱节点：`M1-02`
> 前置：[面力与应力张量](01-traction-and-stress-tensor.md)
> 后续：[主应力](03-principal-stress.md)

---

## 0. 结论先行

对二维应力状态 $(\sigma_x,\sigma_y,\tau_{xy})$，法线从 $x$ 轴逆时针量起角度为 $\theta$ 的斜截面上：

$$
\boxed{\sigma_n = \frac{\sigma_x+\sigma_y}{2}+\frac{\sigma_x-\sigma_y}{2}\cos 2\theta+\tau_{xy}\sin 2\theta}
$$

$$
\boxed{\tau = -\frac{\sigma_x-\sigma_y}{2}\sin 2\theta+\tau_{xy}\cos 2\theta}
$$

三点推论：

1. **正应力 + 剪应力不是把原载荷"分掉"**。它们的平方和等于面力大小的平方： $\sigma_n^2+\tau^2=|\mathbf{t}|^2$。
   投影不改变面力，只是把同一个矢量换了一个分解方式。
2. **角度会翻倍**（ $2\theta$）。这是后面莫尔圆能画成圆的根源。
3. 单向拉伸 100 MPa 时， $\theta=45^\circ$ 的截面上正应力 50 MPa、剪应力 50 MPa。
   这不是巧合： $\cos 90^\circ=0$，一半项消失，剩下 $100/2$。

---

## 1. 前提与符号

- 二维应力状态（平面应力或平面应变，两者在此题中公式相同，见 §4 误解 3）。
- $\theta$：**截面法线**与 $x$ 轴夹角，逆时针为正。
- 取 $\mathbf{n}=(\cos\theta,\ \sin\theta)$，面内与 $\mathbf{n}$ 垂直的方向取 $\mathbf{m}=(-\sin\theta,\ \cos\theta)$。

> 注意：很多教材的 $\theta$ 是从 $x$ 轴量到**截面本身**（与法线差 $90^\circ$），
> 会导致公式中 $\sigma_x$ 与 $\sigma_y$ 互换。读别处的书时先看它的图。

---

## 2. 推导

### 2.1 方法一：张量投影（推荐）

**步骤 1｜写出面力矢量。**
由 $t_i=\sigma_{ij}n_j$，二维展开：

$$
t_1 = \sigma_x\cos\theta+\tau_{xy}\sin\theta
$$

$$
t_2 = \tau_{xy}\cos\theta+\sigma_y\sin\theta
$$

**依据**：Cauchy 公式 + 对称性 $\tau_{xy}=\tau_{yx}$。

**步骤 2｜正应力 = 面力在法向上的投影。**

$$
\sigma_n = \mathbf{t}\cdot\mathbf{n} = t_1\cos\theta+t_2\sin\theta
$$

代入：

$$
\sigma_n=\sigma_x\cos^2\theta+2\tau_{xy}\sin\theta\cos\theta+\sigma_y\sin^2\theta
$$

**依据**：投影的坐标表达式；两个 $\tau$ 项合并（ $\tau_{xy}=\tau_{yx}$）。

**步骤 3｜用倍角公式化简。**
代入 $\cos^2\theta=\frac{1+\cos2\theta}{2}$， $\sin^2\theta=\frac{1-\cos2\theta}{2}$， $2\sin\theta\cos\theta=\sin2\theta$：

$$
\sigma_n=\frac{\sigma_x+\sigma_y}{2}+\frac{\sigma_x-\sigma_y}{2}\cos2\theta+\tau_{xy}\sin2\theta
$$

**依据**：三角恒等式，纯代数。

**步骤 4｜剪应力 = 面力在面内的投影。**

$$
\tau=\mathbf{t}\cdot\mathbf{m}= t_1(-\sin\theta)+t_2\cos\theta
$$

代入 $t_1,t_2$：

$$
\tau=-\sigma_x\sin\theta\cos\theta-\tau_{xy}\sin^2\theta+\tau_{xy}\cos^2\theta+\sigma_y\sin\theta\cos\theta
$$

整理：

$$
\tau=-\frac{\sigma_x-\sigma_y}{2}\sin2\theta+\tau_{xy}\cos2\theta
$$

**依据**：同上，用 $\cos^2\theta-\sin^2\theta=\cos 2\theta$。

> **符号约定**：本式给出的 $\tau$ 以使微元逆时针转为正。
> 若只关心大小，取绝对值即可；莫尔圆作图中习惯把剪应力画在纵轴正方向，需注意符号约定。

### 2.2 方法二：楔形微元静力平衡（独立验证）

取一个直角边分别平行于 $x,y$ 轴、斜边法向为 $\mathbf{n}$ 的三角形微元，
斜边面积 $A$，厚度 $t$（平面问题取单位厚度）。

**步骤 1｜三个面的面积。**
斜边在 $y$ 面上的投影为 $A\cos\theta$，在 $x$ 面上的投影为 $A\sin\theta$。

**步骤 2｜沿法向列力平衡。**
斜边上面力分量 $\sigma_n,\tau$；两个直角边上分别有 $\sigma_x$（沿 $x$）与 $\sigma_y$（沿 $y$），以及剪应力 $\tau_{xy}$。

沿法向投影：

$$
\sigma_n A-\sigma_x (A\cos\theta)\cos\theta-\sigma_y(A\sin\theta)\sin\theta
-2\tau_{xy}(A\sin\theta)\cos\theta = 0
$$

**步骤 3｜约去 $A$，整理。**

$$
\sigma_n=\sigma_x\cos^2\theta+\sigma_y\sin^2\theta+2\tau_{xy}\sin\theta\cos\theta
$$

与方法一**完全一致**。✓

> 两条独立路径给出同一结果，才说明推导没有算错。这是本库"数值/独立校验"要求的来源。

---

## 3. 数值检验

### 3.1 单轴拉伸 100 MPa， $\theta=45^\circ$

$\sigma_x=100$， $\sigma_y=0$， $\tau_{xy}=0$：

$$
\sigma_n=50+50\cos90^\circ = 50\ \mathrm{MPa}
$$

$$
\tau=-50\sin90^\circ=-50\ \mathrm{MPa},\qquad |\tau|=50\ \mathrm{MPa}
$$

**核对面力大小**：

$$
t_1=100\cos45^\circ=70.71,\quad t_2=0,\quad |\mathbf{t}|=70.71
$$

$$
\sqrt{\sigma_n^2+\tau^2}=\sqrt{50^2+50^2}=70.71\ \checkmark
$$

**几何解释**：
面力矢量始终沿 $x$ 方向、大小 70.71（因为它只是被"投影"到斜面上，分量 $t_1$ 被 $\cos\theta$ 缩放）。
把它拆成"垂直于斜面"和"沿斜面"两个方向，正好各 50。

> 图 2：左图为 45° 斜截面上的面力分解；右图为对应的莫尔圆位置（ $\theta=45^\circ$ 对应圆心角 $90^\circ$）。
> 文件：`assets/figures/fig02-oblique-plane-45.png`

### 3.2 特殊角度自检表

| $\theta$ | 截面 | $\sigma_n$ | $\tau$ |
|---|---|---|---|
| $0^\circ$ | 垂直 $x$ 轴 | $\sigma_x=100$ | $0$ |
| $30^\circ$ | 斜 30° | $75$ | $-43.3$ |
| $45^\circ$ | 斜 45° | $50$ | $-50$ |
| $60^\circ$ | 斜 60° | $25$ | $-43.3$ |
| $90^\circ$ | 垂直 $y$ 轴（平行棒轴） | $0$ | $0$ |

$\theta=90^\circ$ 时 $\sigma_n=0$、 $\tau=0$：这正是"棒侧面是自由表面"的数学表达。✓

### 3.3 纯剪状态

$\sigma_x=\sigma_y=0$， $\tau_{xy}=50$：

$$
\sigma_n=50\sin2\theta,\qquad \tau=50\cos2\theta
$$

$\theta=0$： $\sigma_n=0,\tau=50$； $\theta=45^\circ$： $\sigma_n=50,\tau=0$。
**纯剪的 45° 截面上没有剪应力**——这就是主方向。见下一篇。

---

## 4. 常见误解

**误解 1："100 MPa 被分成了 50 + 50 两份力。"**

为什么错： $50^2+50^2=5000 \lt  100^2$。分量平方和小于合矢量平方，说明这是**同一矢量的两次投影**，不是"分成两份"。
面力的法向分量与切向分量**不是**两个独立的力，而是同一个面力在两组正交方向上的坐标。

正确说法：载荷没有变，变的是截面；50/50 是"70.71 MPa 面力"在两个方向上的投影。

**误解 2："正应力一定小于原来的拉应力。"**

为什么错： $\sigma_n$ 在 $\theta=0$ 处等于 100，在 $45^\circ$ 处等于 50，在 $90^\circ$ 处为 0。
若存在双向应力， $\sigma_n$ 可能大于任一单向分量（例如 $\sigma_x=\sigma_y=100$ 时为 100，任何方向都是 100）。

**误解 3："平面应力与平面应变公式不同。"**

为什么错：上述公式只依赖"某个方向上的应力分量为零"，与位移约束无关。
平面应力（ $\sigma_z=0$）与平面应变（ $\varepsilon_z=0$，但 $\sigma_z\ne0$）在 $x$– $y$ 面内的**应力变换公式相同**，
区别只在 $\sigma_z$ 的取值和本构关系。

---

## 5. 自测

1. $\sigma_x=50,\sigma_y=50,\tau_{xy}=0$。求任意 $\theta$ 下的 $\sigma_n,\tau$。
2. $\sigma_x=80,\sigma_y=20,\tau_{xy}=40$。求 $\theta=30^\circ$ 截面上的 $\sigma_n,\tau$。
3. 证明：若 $|\tau_{xy}|$ 足够大，即使 $\sigma_x=\sigma_y=0$，也存在一个截面上 $\sigma_n=0$ 而 $\tau\ne0$。
4. 一个学生把 $\theta$ 量到"截面本身"而不是"法线"，他的公式会差多少？请具体说明。
5. 用 3.2 的表解释：为什么单向拉伸时，平行于棒轴的纵向截面"什么都不承受"？

<details>
<summary>参考答案</summary>

1. $\sigma_n=50$， $\tau=0$。各向同性面内等双轴拉伸：任何方向都是主方向，无剪应力。
2. $\sigma_n=50+30\cos60^\circ+40\sin60^\circ=50+15+34.64=99.64$ MPa；
   $\tau=-30\sin60^\circ+40\cos60^\circ=-25.98+20=-5.98$ MPa。
3. 取 $\sigma_x=\sigma_y=0$： $\sigma_n=\tau_{xy}\sin2\theta$。令 $\theta=0$，则 $\sigma_n=0$ 而 $\tau=\tau_{xy}\ne0$。
   这就是纯剪状态在原坐标面上的表现。
4. 差 $90^\circ$，即 $\theta\to\theta+90^\circ$： $\cos2\theta$、 $\sin2\theta$ 均反号，
   相当于把 $\sigma_x$ 与 $\sigma_y$ 互换、 $\tau$ 反号。数值会完全不同（例如 45° 时 50 vs 0）。
5. $\theta=90^\circ$： $\sigma_n=\sigma_x\cos^2 90^\circ=0$， $\tau=-\tfrac{\sigma_x}{2}\sin180^\circ=0$。
   纵向截面是自由表面，没有材料跨越它传递内力。
</details>

---

下一篇 → [主应力：为什么切应力能消失](03-principal-stress.md)
