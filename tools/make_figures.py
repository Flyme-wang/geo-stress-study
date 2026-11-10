#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
geo-stress-study 教学图件生成脚本

运行：
    python tools/make_figures.py

输出：assets/figures/*.png

约定：
- 拉为正；正应力用绿色，剪应力用红色，面力/参考线用蓝色。
- 所有数值均为教科书算例（单轴 100 MPa / 纯剪 50 MPa），便于与正文交叉核对。
- 图件请勿手改；修改本脚本后重跑。
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Wedge, FancyArrowPatch, Rectangle, Arc

# ---------------------------------------------------------------- 全局样式
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "assets", "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "font.size": 11,
    "axes.linewidth": 0.9,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

GREEN = "#2E7D32"   # 正应力
RED = "#C62828"     # 剪应力
BLUE = "#1565C0"    # 面力 / 参考线
GRAY = "#5F6B73"

C30 = np.cos(np.radians(30))
S30 = np.sin(np.radians(30))


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path)
    plt.close(fig)
    print("saved", os.path.normpath(path))


def arrow(ax, p0, p1, color, lw=2.0, ms=13, ls="-", zorder=6):
    ax.add_patch(FancyArrowPatch(
        np.asarray(p0, float), np.asarray(p1, float),
        arrowstyle="-|>", mutation_scale=ms, color=color, lw=lw,
        linestyle=ls, zorder=zorder, shrinkA=0, shrinkB=0))


# ============================================================ 图 1 应力微元
def fig01():
    """等轴测投影：u=(y-x)cos30, v=z-(x+y)sin30"""
    def P(x, y, z):
        return np.array([(y - x) * C30, z - (x + y) * S30])

    ex = P(1, 0, 0) - P(0, 0, 0)   # +x 屏方向
    ey = P(0, 1, 0) - P(0, 0, 0)
    ez = np.array([0.0, 1.0])

    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    ax.set_aspect("equal")
    ax.axis("off")

    faces = {
        "top":   ([(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)], "#D8EAFB"),
        "left":  ([(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)], "#A9CCEF"),
        "right": ([(1, 1, 0), (0, 1, 0), (0, 1, 1), (1, 1, 1)], "#8CB6DE"),
    }
    for vs, c in faces.values():
        ax.add_patch(Polygon([P(*v) for v in vs], closed=True, facecolor=c,
                             edgecolor="#37474F", lw=1.1, zorder=2))

    # ---- 正应力（绿色）
    c_left = P(1, 0.5, 0.5)
    arrow(ax, c_left, c_left + 0.55 * ex, GREEN, lw=2.4)
    ax.annotate(r"$\sigma_{11}$", c_left + 0.55 * ex,
                xytext=(-2, -13), textcoords="offset points",
                color=GREEN, fontsize=13, ha="center", va="center")

    c_right = P(0.5, 1, 0.5)
    arrow(ax, c_right, c_right + 0.55 * ey, GREEN, lw=2.4)
    ax.annotate(r"$\sigma_{22}$", c_right + 0.55 * ey,
                xytext=(10, -8), textcoords="offset points",
                color=GREEN, fontsize=13, ha="center", va="center")

    c_top = P(0.5, 0.5, 1)
    arrow(ax, c_top, c_top + 0.60 * ez, GREEN, lw=2.4)
    ax.annotate(r"$\sigma_{33}$", c_top + 0.60 * ez,
                xytext=(10, 0), textcoords="offset points",
                color=GREEN, fontsize=13, ha="left", va="center")

    # ---- 剪应力（红色）
    b = P(1, 0.30, 0.22)
    arrow(ax, b, b + 0.30 * ey, RED, lw=2.0, ms=11)
    ax.annotate(r"$\tau_{12}$", b + 0.30 * ey, xytext=(2, -13),
                textcoords="offset points", color=RED, fontsize=12, ha="center")
    arrow(ax, b, b + 0.28 * ez, RED, lw=2.0, ms=11)
    ax.annotate(r"$\tau_{13}$", b + 0.28 * ez, xytext=(-10, 2),
                textcoords="offset points", color=RED, fontsize=12, ha="right")

    b2 = P(0.30, 1, 0.22)
    arrow(ax, b2, b2 + 0.30 * ex, RED, lw=2.0, ms=11)
    ax.annotate(r"$\tau_{21}$", b2 + 0.30 * ex, xytext=(-2, -13),
                textcoords="offset points", color=RED, fontsize=12, ha="center")
    arrow(ax, b2, b2 + 0.28 * ez, RED, lw=2.0, ms=11)
    ax.annotate(r"$\tau_{23}$", b2 + 0.28 * ez, xytext=(8, 2),
                textcoords="offset points", color=RED, fontsize=12, ha="left")

    b3 = P(0.62, 0.30, 1)
    arrow(ax, b3, b3 + 0.26 * ex, RED, lw=2.0, ms=11)
    ax.annotate(r"$\tau_{31}$", b3 + 0.26 * ex, xytext=(-6, 6),
                textcoords="offset points", color=RED, fontsize=12, ha="right")
    arrow(ax, b3, b3 + 0.26 * ey, RED, lw=2.0, ms=11)
    ax.annotate(r"$\tau_{32}$", b3 + 0.26 * ey, xytext=(4, 8),
                textcoords="offset points", color=RED, fontsize=12, ha="left")

    # ---- 坐标三面角（画在左上角，避免与立方体重叠）
    o = np.array([-1.55, 0.35])
    for d, lab, off in [(ex, r"$x_1$", (-8, -8)),
                        (ey, r"$x_2$", (4, -8)),
                        (ez, r"$x_3$", (2, 4))]:
        tip = o + 0.30 * d
        ax.plot([o[0], tip[0]], [o[1], tip[1]], color="#37474F", lw=1.0)
        ax.annotate("", xy=tip, xytext=o, arrowprops=dict(
            arrowstyle="-|>", color="#37474F", lw=1.0, mutation_scale=9))
        ax.annotate(lab, tip, xytext=off, textcoords="offset points",
                    fontsize=11, color="#37474F")

    # ---- 图例
    ax.plot([], [], color=GREEN, lw=2.4, label="正应力 $\\sigma_{ij}$（沿面法向）")
    ax.plot([], [], color=RED, lw=2.0, label="剪应力 $\\tau_{ij}$（面内方向）")
    ax.legend(loc="upper right", fontsize=9.5, frameon=True, framealpha=0.95,
              borderpad=0.6, handlelength=1.8)

    ax.set_title("三维微元三个可见面上的应力分量（拉为正）", fontsize=13, pad=8)
    ax.text(0.0, -0.06,
            "下标含义：$\\sigma_{ij}$ 的第 1 个下标为该面法向，第 2 个为力的方向；"
            "对称性给出 $\\tau_{12}=\\tau_{21}$",
            transform=ax.transAxes, fontsize=9.5, color=GRAY, ha="center")
    ax.set_xlim(-2.0, 1.8)
    ax.set_ylim(-1.35, 1.45)
    save(fig, "fig01-stress-cube.png")


# ==================================================== 图 2 45° 斜截面
def fig02():
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.6))
    ax = axes[0]
    ax.set_aspect("equal")
    ax.axis("off")

    # 三角微元：直角边沿 x、y，斜边为 45° 截面
    tri = Polygon([[0, 0], [1.15, 0], [0, 1.15]], closed=True,
                  facecolor="#E8F1FB", edgecolor="#37474F", lw=1.3, zorder=2)
    ax.add_patch(tri)

    c = np.array([1.15, 1.15]) / 2.0      # 斜边中点
    n = np.array([1, 1]) / np.sqrt(2)     # 法向
    t_dir = np.array([1, -1]) / np.sqrt(2)

    arrow(ax, c, c + 0.62 * np.array([1, 0]), BLUE, lw=2.8, ms=15)
    ax.annotate(r"$\mathbf{t}=70.7$ MPa", c + np.array([0.30, 0.02]),
                color=BLUE, fontsize=11, va="bottom")

    arrow(ax, c, c + 0.46 * n, GREEN, lw=2.4)
    ax.annotate(r"$\sigma_n=50$", c + 0.50 * n, xytext=(6, 4),
                textcoords="offset points", color=GREEN, fontsize=11,
                ha="left", va="bottom")

    arrow(ax, c, c + 0.46 * t_dir, RED, lw=2.4)
    ax.annotate(r"$\tau=50$", c + 0.50 * t_dir, xytext=(6, -4),
                textcoords="offset points", color=RED, fontsize=11,
                ha="left", va="top")

    # 直角标记
    d = 0.07
    ax.plot([c[0] + d * n[0], c[0] + d * (n[0] + t_dir[0]),
             c[0] + d * t_dir[0]],
            [c[1] + d * n[1], c[1] + d * (n[1] + t_dir[1]),
             c[1] + d * t_dir[1]], color=GRAY, lw=0.9)

    # 45° 角标（在左上顶点，介于竖边与斜边之间）
    ax.add_patch(Arc((0, 1.15), 0.42, 0.42, theta1=-90, theta2=-45,
                     color=GRAY, lw=0.9))
    ax.annotate(r"$45^\circ$", (0.13, 0.88), fontsize=10, color=GRAY)
    ax.annotate("45° 斜截面", (0.10, 0.60), fontsize=10, color=GRAY,
                ha="left", va="bottom")

    ax.set_xlim(-0.35, 2.15)
    ax.set_ylim(-0.35, 1.75)
    ax.set_title("单轴 100 MPa：45° 斜截面上的面力分解", fontsize=12)

    # ---- 右：莫尔圆
    ax = axes[1]
    C, R = 50.0, 50.0
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(C + R * np.cos(th), R * np.sin(th), color="#37474F", lw=1.6)
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.axvline(0, color=GRAY, lw=0.8)

    p0 = (100.0, 0.0)
    p1 = (50.0, -50.0)
    p2 = (0.0, 0.0)
    for (x, y), lab, off in [(p0, r"$\theta=0^\circ$", (10, 8)),
                             (p1, r"$\theta=45^\circ$", (12, -14)),
                             (p2, r"$\theta=90^\circ$", (-64, 6))]:
        ax.plot(x, y, "o", color=RED, ms=7, zorder=6)
        ax.annotate(lab, (x, y), xytext=off, textcoords="offset points",
                    fontsize=10, color="#37474F")

    ax.add_patch(Wedge((C, 0), 24, -90, 0, facecolor="#EF6C00",
                       alpha=0.18, edgecolor="none"))
    ax.add_patch(Arc((C, 0), 48, 48, theta1=-90, theta2=0,
                     color="#EF6C00", lw=1.4))
    ax.plot([C, p0[0]], [0, 0], color="#EF6C00", lw=0.9, ls=":")
    ax.plot([C, p1[0]], [0, p1[1]], color="#EF6C00", lw=0.9, ls=":")
    ax.annotate(r"圆上转 $2\theta=90^\circ$", (C + 2, -46),
                color="#EF6C00", fontsize=10)

    ax.plot([C], [0], "o", color=GRAY, ms=4)
    ax.annotate("圆心 $C=50$", (C, 0), xytext=(10, 10),
                textcoords="offset points", fontsize=10, color=GRAY)
    ax.set_xlim(-28, 128)
    ax.set_ylim(-72, 66)
    ax.set_xlabel(r"$\sigma_n$ / MPa")
    ax.set_ylabel(r"$\tau$ / MPa")
    ax.set_aspect("equal")
    ax.set_title("对应的莫尔圆：实物转 45°，圆上转 90°", fontsize=12)
    ax.grid(alpha=0.22)

    fig.tight_layout()
    save(fig, "fig02-oblique-plane-45.png")


# ==================================================== 图 3 三维莫尔圆
def fig03():
    s1, s2, s3 = 100.0, 40.0, -20.0
    fig, ax = plt.subplots(figsize=(8.8, 6.4))
    th = np.linspace(0, 2 * np.pi, 600)

    def circ(a, b, **kw):
        C, R = (a + b) / 2, (b - a) / 2
        ax.plot(C + R * np.cos(th), R * np.sin(th), **kw)

    # 容许区：大圆内、两小圆外
    Cb, Rb = (s1 + s3) / 2, (s1 - s3) / 2
    ax.fill(Cb + Rb * np.cos(th), Rb * np.sin(th), color="#DCEBFA", zorder=0)
    for a, b in [(s1, s2), (s2, s3)]:
        C, R = (a + b) / 2, (b - a) / 2
        ax.fill(C + R * np.cos(th), R * np.sin(th), color="white", zorder=1)

    circ(s1, s3, color="#263238", lw=1.9, zorder=3)
    circ(s1, s2, color="#546E7A", lw=1.3, ls="--", zorder=3)
    circ(s2, s3, color="#546E7A", lw=1.3, ls="--", zorder=3)

    ax.axhline(0, color=GRAY, lw=0.9, zorder=4)
    box = dict(facecolor="white", edgecolor="none", pad=1.5)
    for s, lab, dy in [(s1, r"$\sigma_1=100$", -30),
                       (s2, r"$\sigma_2=40$", -30),
                       (s3, r"$\sigma_3=-20$", -30)]:
        ax.plot([s], [0], "o", color="#263238", ms=5, zorder=6)
        ax.annotate(lab, (s, 0), xytext=(0, dy), textcoords="offset points",
                    ha="center", fontsize=10, bbox=box, zorder=7)

    Cmax, Rmax = (s1 + s3) / 2, (s1 - s3) / 2
    ax.plot([Cmax], [Rmax], "o", color=RED, ms=8, zorder=6)
    ax.plot([Cmax, Cmax], [-Rmax, Rmax], ls=":", color=RED, lw=1.1, zorder=2)
    ax.annotate(r"$\tau_{\max}=(\sigma_1-\sigma_3)/2=60$ MPa",
                (Cmax, Rmax), xytext=(12, -4), textcoords="offset points",
                color=RED, fontsize=11, va="center")
    ax.annotate(r"$\sigma_n=40$ MPa", (Cmax, 18), xytext=(-8, 0),
                textcoords="offset points", color=RED, fontsize=10, ha="right")

    ax.set_xlabel(r"$\sigma_n$ / MPa")
    ax.set_ylabel(r"$\tau$ / MPa")
    ax.set_title("三维应力状态的三个莫尔圆（阴影为可能的截面状态）", fontsize=12)
    ax.set_aspect("equal")
    ax.grid(alpha=0.22)
    ax.set_xlim(-42, 138)
    ax.set_ylim(-76, 76)
    ax.text(0.02, 0.03, "两小圆内部不是任何截面的状态；\n截面状态只能落在阴影曲边三角形内",
            transform=ax.transAxes, fontsize=9.5, color=GRAY)
    save(fig, "fig03-mohr-circles-3d.png")


# ==================================================== 图 4 强度准则
def fig04():
    fig, ax = plt.subplots(figsize=(7.8, 7.0))
    a = 1.0

    ax.add_patch(Rectangle((-1, -1), 2, 2, fill=False, edgecolor=BLUE,
                           lw=2.2, label="Rankine（最大主应力）"))

    t = np.linspace(0, 2 * np.pi, 500)
    A, B = np.sqrt(2) * a, np.sqrt(2 / 3) * a
    x = (A * np.cos(t) + B * np.sin(t)) / np.sqrt(2)
    y = (-A * np.cos(t) + B * np.sin(t)) / np.sqrt(2)
    ax.plot(x, y, color="#E65100", lw=2.4, label="von Mises（畸变能）")

    hexa = [(1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1), (1, 0), (1, 1)]
    ax.plot([p[0] for p in hexa], [p[1] for p in hexa], ls="--",
            color="#2E7D32", lw=2.0, label="Tresca（最大剪应力）")

    ax.plot([0, 1], [0, 0], color=GRAY, lw=1.0, ls=":")
    ax.annotate("单轴拉伸", (0.52, 0), xytext=(0, 9), textcoords="offset points",
                fontsize=9, color=GRAY, ha="center")
    ax.annotate("纯剪：两准则相差 15.5%", (0.54, -0.60),
                xytext=(0, -14), textcoords="offset points", fontsize=9,
                color=GRAY, ha="center")
    ax.plot([0.5, 0.577], [-0.5, -0.577], "o", color=GRAY, ms=4, ls="none")

    ax.axhline(0, color=GRAY, lw=0.8)
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_xlabel(r"$\sigma_2/\sigma_y$")
    ax.set_ylabel(r"$\sigma_1/\sigma_y$")
    ax.set_aspect("equal")
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.grid(alpha=0.22)
    ax.set_title("平面应力下的三种屈服准则（$\\sigma_3=0$）", fontsize=12)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3,
              fontsize=9.5, frameon=False)
    ax.text(0.02, 0.03,
            "四象限对称；von Mises 外接 Tresca 六边形，\n二者差异最大处为纯剪（约 15.5%）",
            transform=ax.transAxes, fontsize=9.5, color=GRAY)
    save(fig, "fig04-yield-criteria.png")


# ==================================================== 图 5 单轴应力路径
def fig05():
    fig, ax = plt.subplots(figsize=(7.8, 6.2))
    C, R = 50.0, 50.0
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(C + R * np.cos(th), R * np.sin(th), color="#263238", lw=1.6)
    ax.fill(C + R * np.cos(th), R * np.sin(th), color="#E8F1FB", zorder=0)

    angles = [0, 15, 30, 45, 60, 75, 90]
    xs, ys = [], []
    for a_ in angles:
        x = C + R * np.cos(np.radians(2 * a_))
        y = -R * np.sin(np.radians(2 * a_))
        xs.append(x)
        ys.append(y)
        ax.plot(x, y, "o", ms=5, color="#546E7A", zorder=5)

    for a_, lab, off in [(0, r"$\theta=0^\circ$", (10, 6)),
                         (45, r"$\theta=45^\circ$", (0, -22)),
                         (90, r"$\theta=90^\circ$", (10, 12))]:
        x = C + R * np.cos(np.radians(2 * a_))
        y = -R * np.sin(np.radians(2 * a_))
        ax.plot(x, y, "o", ms=7, color=RED, zorder=6)
        ax.annotate(lab, (x, y), xytext=off, textcoords="offset points",
                    fontsize=10, color="#263238")

    ax.plot([C, C + R], [0, 0], color=GRAY, lw=0.9, ls=":")
    ax.plot([C, C], [0, -R], color=GRAY, lw=0.9, ls=":")
    ax.add_patch(Arc((C, 0), 44, 44, theta1=-90, theta2=0,
                     color=GRAY, lw=1.1, ls=":"))
    ax.annotate(r"$2\theta=90^\circ$", (C + 24, -16), color=GRAY, fontsize=10)

    ax.axhline(0, color=GRAY, lw=0.9)
    ax.axvline(0, color=GRAY, lw=0.9)
    ax.set_xlim(-28, 126)
    ax.set_ylim(-70, 62)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\sigma_n$ / MPa")
    ax.set_ylabel(r"$\tau$ / MPa")
    ax.set_title("单轴拉伸 100 MPa：截面转过 90°，圆上转过 180°", fontsize=12)
    ax.grid(alpha=0.22)
    ax.text(0.03, 0.04,
            "截面角 $\\theta$ 从 $0^\\circ\\to45^\\circ\\to90^\\circ$，"
            "对应点从 (100,0) 经 (50,-50) 到 (0,0)",
            transform=ax.transAxes, fontsize=9.5, color=GRAY)
    save(fig, "fig05-stress-path-uniaxial.png")


# ==================================================== 图 6 断层类型
def fig06():
    """三类断层。第 1、3 幅为剖面，第 2 幅为俯视（走滑需俯视才看得出错动）。"""
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.6))
    xmin, xmax, ytop, ybot = -1.25, 1.25, 0.80, -1.05

    # ---------------------------------------------------- 剖面：正断层
    ax = axes[0]
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-2.10, 1.30)
    ax.add_patch(Rectangle((xmin, ybot), xmax - xmin, ytop - ybot,
                           facecolor="#F4F1EC", edgecolor="none", zorder=1))
    ax.plot([xmin, xmax], [ytop, ytop], color="#37474F", lw=1.6, zorder=3)
    for x in np.arange(xmin + 0.06, xmax, 0.24):
        ax.plot([x, x - 0.10], [ytop, ytop + 0.11], color="#37474F", lw=0.8,
                zorder=3)

    th = np.radians(60)
    ts = np.linspace(-3, 3, 1500)
    X, Y = np.cos(th) * ts, -np.sin(th) * ts
    ok = (X >= xmin) & (X <= xmax) & (Y >= ybot) & (Y <= ytop)
    ax.plot(X[ok], Y[ok], color="#6D4C41", lw=3.2, zorder=4,
            solid_capstyle="butt")
    ux, uy = -np.cos(th), np.sin(th)      # 向断层上倾方向
    arrow(ax, (-0.30, -0.12), (-0.30, -0.12) + 0.34 * np.array([ux, uy]),
          "#8D6E63", lw=1.9, ms=11)
    arrow(ax, (0.30, -0.12), (0.30, -0.12) - 0.34 * np.array([ux, uy]),
          "#8D6E63", lw=1.9, ms=11)
    ax.annotate("上盘下降", (0.62, -0.62), fontsize=9.5, color="#6D4C41")

    o = np.array([0.0, -1.55])
    arrow(ax, o + [0, -0.40], o + [0, 0.32], RED, lw=2.2, ms=12)
    ax.annotate(r"$\sigma_1$（竖直）", o + [0, 0.36], xytext=(6, 0),
                textcoords="offset points", color=RED, fontsize=10.5,
                va="bottom")
    arrow(ax, o + [0.02, -0.05], o + [-0.62, -0.05], GREEN, lw=2.0, ms=12)
    ax.annotate(r"$\sigma_3$", o + [-0.66, -0.05], xytext=(-6, -4),
                textcoords="offset points", color=GREEN, fontsize=10.5,
                ha="right", va="center")
    ax.set_title("正断层（伸展）", fontsize=12)

    # ---------------------------------------------------- 俯视：走滑断层
    ax = axes[1]
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-2.10, 1.30)
    ax.add_patch(Rectangle((xmin, ybot), xmax - xmin, ytop - ybot,
                           facecolor="#F4F1EC", edgecolor="#B0A89C",
                           lw=0.8, zorder=1))
    ax.annotate("俯视（从上方看）", (0, ytop + 0.18), ha="center",
                fontsize=9.5, color=GRAY)
    ax.plot([0, 0], [ybot, ytop], color="#6D4C41", lw=3.2, zorder=4,
            solid_capstyle="butt")
    arrow(ax, (-0.42, -0.25), (-0.42, 0.42), "#8D6E63", lw=2.0, ms=12)
    arrow(ax, (0.42, 0.25), (0.42, -0.42), "#8D6E63", lw=2.0, ms=12)
    ax.annotate("左行错动", (-0.42, 0.52), ha="center", fontsize=9.5,
                color="#6D4C41")

    o = np.array([0.0, -1.55])
    arrow(ax, o + [-0.10, -0.05], o + [-0.62, -0.05], RED, lw=2.2, ms=12)
    ax.annotate(r"$\sigma_1$", o + [-0.68, -0.05], xytext=(-6, -4),
                textcoords="offset points", color=RED, fontsize=10.5,
                ha="right", va="center")
    arrow(ax, o + [0.10, -0.05], o + [0.62, -0.05], GREEN, lw=2.0, ms=12)
    ax.annotate(r"$\sigma_3$", o + [0.68, -0.05], xytext=(6, -4),
                textcoords="offset points", color=GREEN, fontsize=10.5,
                va="center")
    ax.annotate(r"$\sigma_2$ 竖直（指向纸外）", (0, -1.92), ha="center",
                fontsize=9.5, color=GRAY)
    ax.set_title("走滑断层", fontsize=12)

    # ---------------------------------------------------- 剖面：逆冲断层
    ax = axes[2]
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-2.10, 1.30)
    ax.add_patch(Rectangle((xmin, ybot), xmax - xmin, ytop - ybot,
                           facecolor="#F4F1EC", edgecolor="none", zorder=1))
    ax.plot([xmin, xmax], [ytop, ytop], color="#37474F", lw=1.6, zorder=3)
    for x in np.arange(xmin + 0.06, xmax, 0.24):
        ax.plot([x, x - 0.10], [ytop, ytop + 0.11], color="#37474F", lw=0.8,
                zorder=3)

    th = np.radians(30)
    x0 = 0.55
    ts = np.linspace(-3, 3, 1500)
    X, Y = x0 + np.cos(th) * ts, -np.sin(th) * ts
    ok = (X >= xmin) & (X <= xmax) & (Y >= ybot) & (Y <= ytop)
    ax.plot(X[ok], Y[ok], color="#6D4C41", lw=3.2, zorder=4,
            solid_capstyle="butt")
    ux, uy = -np.cos(th), np.sin(th)
    arrow(ax, (-0.55, -0.35), (-0.55, -0.35) + 0.36 * np.array([ux, uy]),
          "#8D6E63", lw=1.9, ms=11)
    arrow(ax, (0.05, -0.35), (0.05, -0.35) - 0.36 * np.array([ux, uy]),
          "#8D6E63", lw=1.9, ms=11)
    ax.annotate("上盘抬升", (-0.55, 0.30), fontsize=9.5, color="#6D4C41",
                ha="center")

    o = np.array([0.0, -1.55])
    arrow(ax, o + [-0.10, -0.05], o + [-0.62, -0.05], RED, lw=2.2, ms=12)
    ax.annotate(r"$\sigma_1$", o + [-0.68, -0.05], xytext=(-6, -4),
                textcoords="offset points", color=RED, fontsize=10.5,
                ha="right", va="center")
    arrow(ax, o + [0.45, 0.30], o + [0.45, -0.48], GREEN, lw=2.0, ms=12)
    ax.annotate(r"$\sigma_3$（竖直）", o + [0.45, -0.52], xytext=(4, -4),
                textcoords="offset points", color=GREEN, fontsize=10.5,
                va="top")
    ax.set_title("逆冲断层（挤压）", fontsize=12)

    fig.suptitle("Anderson 断层分类：竖直方向对应哪个主应力", fontsize=13, y=1.02)
    fig.tight_layout()
    save(fig, "fig06-anderson-faults.png")


if __name__ == "__main__":
    fig01()
    fig02()
    fig03()
    fig04()
    fig05()
    fig06()
