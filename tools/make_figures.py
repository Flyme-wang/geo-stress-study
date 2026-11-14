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
from matplotlib.patches import (Polygon, Wedge, FancyArrowPatch, Rectangle, Arc,
                                FancyBboxPatch, Ellipse)
import element_lib as el

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


# ======================================== 图 7 单单元算例：模型、位移与应力
def _q4_square(E, nu, t, F):
    """1×1 正方形 Q4（节点 1-2-3-4 逆时针），四个高斯点积分。

    固定节点 1、2 的全部自由度；返回 8 个自由度的位移、以及可在
    单元内任意点求值的场函数 (ε, σ)。
    """
    coords = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    D = E / (1 - nu ** 2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])
    g = 1 / np.sqrt(3)

    def Bmat(xi, eta):
        dN = np.array([
            [-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],
            [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]]) / 4.0
        J = dN @ coords
        dNdX = np.linalg.solve(J, dN)
        B = np.zeros((3, 8))
        B[0, 0::2] = dNdX[0]
        B[1, 1::2] = dNdX[1]
        B[2, 0::2] = dNdX[1]
        B[2, 1::2] = dNdX[0]
        return B, np.linalg.det(J)

    K = np.zeros((8, 8))
    for xi in (-g, g):
        for eta in (-g, g):
            B, detJ = Bmat(xi, eta)
            K += B.T @ D @ B * detJ * t
    u = np.zeros(8)
    u[4:] = np.linalg.solve(K[4:, 4:], F[4:])

    def field(xi, eta):
        B, _ = Bmat(xi, eta)
        return B @ u, D @ B @ u

    return u, field


def fig07():
    import numpy as np
    E, nu, t = 300e9, 0.3, 1.0
    F = np.array([0, 0, 0, 0, 0, -100.0, 0, -100.0])
    u, field = _q4_square(E, nu, t, F)

    xi = np.linspace(-1, 1, 81)
    XI, ETA = np.meshgrid(xi, xi)
    X = (1 + XI) / 2
    Y = (1 + ETA) / 2
    N = np.stack([(1 - XI) * (1 - ETA), (1 + XI) * (1 - ETA),
                  (1 + XI) * (1 + ETA), (1 - XI) * (1 + ETA)], axis=0) / 4

    N0, N1, N2, N3 = N[0], N[1], N[2], N[3]
    ux = (N0 * u[0] + N1 * u[2] + N2 * u[4] + N3 * u[6])
    uy = (N0 * u[1] + N1 * u[3] + N2 * u[5] + N3 * u[7])
    umag = np.hypot(ux, uy)

    sig = np.zeros_like(X)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            sig[i, j] = field(XI[i, j], ETA[i, j])[1][1]

    fig, axes = plt.subplots(1, 4, figsize=(14.4, 4.05))

    # ---- ① 模型
    ax = axes[0]
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-0.30, 1.30); ax.set_ylim(-0.42, 1.60)
    ax.add_patch(Rectangle((0, 0), 1, 1, facecolor="#F4F7FA",
                           edgecolor="#37474F", lw=1.6, zorder=2))
    nodes = {1: (0, 0), 2: (1, 0), 3: (1, 1), 4: (0, 1)}
    off = {1: (-16, -12), 2: (10, -12), 3: (10, 8), 4: (-24, 8)}
    for k, (x, y) in nodes.items():
        ax.plot([x], [y], "o", ms=6.5, mfc="white", mec="#37474F", mew=1.5, zorder=4)
        ax.annotate(str(k), (x, y), xytext=off[k], textcoords="offset points",
                    fontsize=11.5, color="#22303A", zorder=5)
    for x in (0.0, 1.0):
        ax.add_patch(Polygon([[x - 0.10, -0.10], [x + 0.10, -0.10],
                              [x, 0.0]], closed=True,
                             facecolor="#BFC8CE", edgecolor="#37474F", lw=1.0, zorder=3))
        ax.plot([x - 0.13, x + 0.13], [-0.10, -0.10], color="#37474F", lw=1.3, zorder=3)
    for x in (0.0, 1.0):
        arrow(ax, (x, 1.40), (x, 1.03), BLUE, lw=2.2, ms=13)
        ax.annotate("100 N", (x, 1.42), xytext=(0, 2), textcoords="offset points",
                    ha="center", va="bottom", fontsize=10.5, color=BLUE)
    ax.annotate("", xy=(0, -0.26), xytext=(1, -0.26),
                arrowprops=dict(arrowstyle="<->", color=GRAY, lw=1.0))
    ax.text(0.5, -0.36, "1 m", ha="center", va="top", fontsize=10.5, color=GRAY)
    ax.set_title("① 模型：1 个 Q4，固定 1、2", fontsize=11.5, pad=8)
    ax.text(0.5, 0.5, "E = 300 GPa\nν = 0.3，t = 1 m\n平面应力",
            ha="center", va="center", fontsize=10, color="#4A5A66",
            linespacing=1.6, zorder=5)

    # ---- ② ③ 位移云图
    for ax, val, title, unit_note in (
            (axes[1], ux * 1e10, "② $u_x$ 云图", "u_x"),
            (axes[2], umag * 1e10, "③ 总位移 $|u|$ 云图", "|u|")):
        pm = ax.pcolormesh(X, Y, val, shading="gouraud", cmap="jet")
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_xlim(-0.04, 1.04); ax.set_ylim(-0.04, 1.04)
        ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor="#37474F", lw=1.3))
        cb = fig.colorbar(pm, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=8.5)
        cb.set_label(unit_note + r" / $10^{-10}$ m", fontsize=9.5)
        ax.set_title(title, fontsize=11.5, pad=8)

    axes[1].text(0.06, 0.94, "顶边 $\\pm1.338$", fontsize=9.5, color="#22303A",
                 ha="left", va="top", zorder=6,
                 bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.6))
    axes[1].text(0.06, 0.06, "底边 0", fontsize=9.5, color="#22303A",
                 ha="left", va="bottom", zorder=6,
                 bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.6))
    axes[2].text(0.06, 0.94, "角点 6.6051", fontsize=9.5, color="#22303A",
                 ha="left", va="top", zorder=6,
                 bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.6))
    axes[2].text(0.06, 0.06, "底边 0", fontsize=9.5, color="#22303A",
                 ha="left", va="bottom", zorder=6,
                 bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.6))

    # ---- ④ σy 云图
    ax = axes[3]
    pm = ax.pcolormesh(X, Y, sig, shading="gouraud", cmap="jet")
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-0.04, 1.04); ax.set_ylim(-0.04, 1.04)
    ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor="#37474F", lw=1.3))
    cb = fig.colorbar(pm, ax=ax, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=8.5)
    cb.set_label(r"$\sigma_y$ / Pa", fontsize=9.5)
    ax.set_title(r"④ $\sigma_y$ 云图", fontsize=11.5, pad=8)
    ax.text(0.5, 0.52, "中心 $-200$\n$=-200\\,/\\,1$ N/m$^2$",
            ha="center", va="center", fontsize=9.5, color="#22303A", zorder=6,
            bbox=dict(facecolor="white", alpha=0.78, edgecolor="none", pad=2.0))
    ax.text(0.06, 0.06, "底边 $-213$", fontsize=9.5, color="white",
            ha="left", va="bottom", zorder=6)
    ax.text(0.06, 0.94, "顶边 $-187$", fontsize=9.5, color="white",
            ha="left", va="top", zorder=6)

    fig.tight_layout()
    save(fig, "fig07-single-element.png")


# ==================================================== 图 8 有限元程序的四步闭环
def fig08():
    fig, ax = plt.subplots(figsize=(12.6, 4.2))
    ax.set_xlim(0, 12.6); ax.set_ylim(0, 4.2); ax.axis("off")

    steps = [
        ("① 离散", "连续体 → 有限个单元\n节点、单元、自由度编号", "#E8F1FB"),
        ("② 单元刚度", "Ke = ∫ BᵀDB t dA\n形函数 + 雅可比 + 高斯积分", "#EAF4EC"),
        ("③ 组装", "K(I,J) = Σ Ke(i,j)\n按自由度映射“累加”", "#FDF3E3"),
        ("④ 求解", "K u = F\n约束 + 荷载 → 解方程", "#F7E9EA"),
    ]
    w, h = 2.62, 1.72
    for i, (title, body, color) in enumerate(steps):
        x = 0.30 + i * 3.10
        ax.add_patch(FancyBboxPatch(
            (x, 1.48), w, h, boxstyle="round,pad=0.10,rounding_size=0.14",
            facecolor=color, edgecolor="#37474F", lw=1.0, zorder=2))
        ax.text(x + w / 2, 1.48 + h - 0.36, title, ha="center", va="center",
                fontsize=13.5, color="#22303A", zorder=3)
        ax.text(x + w / 2, 1.48 + 0.62, body, ha="center", va="center",
                fontsize=9.8, color="#4A5A66", linespacing=1.7, zorder=3)
        if i < 3:
            arrow(ax, (x + w + 0.08, 2.34), (x + w + 0.38, 2.34), GRAY,
                  lw=1.6, ms=11)

    ax.text(6.3, 0.80, "化整为零 → 局部近似 → 整体组装 → 求解方程",
            ha="center", fontsize=11.5, color=GRAY)
    ax.text(6.3, 0.38, "后处理（应变 ε = Bu，应力 σ = Dε）不改变求解，只改变“看结果”的方式",
            ha="center", fontsize=9.8, color="#7A8792")
    save(fig, "fig08-fem-workflow.png")


# ==================================================== 图 9 节点排序与 det(J)
def fig09():
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.9))
    node_xy = {1: (0, 0), 2: (1, 0), 3: (1, 1), 4: (0, 1)}
    cases = [
        ("正确：1→2→3→4", [1, 2, 3, 4], "det J = +0.25\n（四个高斯点同号）", GREEN),
        ("反向：1→4→3→2", [1, 4, 3, 2], "det J = −0.25\n（四个高斯点同号）", RED),
        ("交叉：1→2→4→3", [1, 2, 4, 3], "det J = −0.25η\n（η = 0 处为 0）", RED),
    ]
    for ax, (title, order, note, color) in zip(axes, cases):
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_xlim(-0.34, 1.34); ax.set_ylim(-0.52, 1.42)
        ax.add_patch(Polygon([node_xy[k] for k in order], closed=True,
                             facecolor="#F4F7FA", edgecolor="#37474F",
                             lw=1.4, zorder=2))
        if "交叉" in title:
            ax.add_patch(Rectangle((0, 0.5), 1, 0.5, facecolor="#F6DEDF",
                                   edgecolor="none", alpha=0.85, zorder=1))
            ax.plot([0, 1], [0.5, 0.5], color=RED, lw=1.0, ls=(0, (4, 3)),
                    zorder=3)
            ax.annotate("η = 0：映射折叠", (0.5, 0.5), xytext=(0, -16),
                        textcoords="offset points", ha="center",
                        fontsize=9, color=RED)
        # 按局部顺序画方向箭头
        for a, b in zip(order, order[1:] + order[:1]):
            p, q = np.array(node_xy[a], float), np.array(node_xy[b], float)
            d = q - p
            arrow(ax, p + 0.10 * d, p + 0.86 * d, "#7A8792", lw=1.2, ms=9)
        for k, (x, y) in node_xy.items():
            local = order.index(k) + 1
            ax.plot([x], [y], "o", ms=9, color=color, zorder=6,
                    markeredgecolor="white", markeredgewidth=1.1)
            dx = -0.10 if x < 0.5 else 0.10
            dy = -0.10 if y < 0.5 else 0.10
            ax.annotate("节点 %d\n（局部 %d）" % (k, local), (x, y),
                        xytext=(dx * 130, dy * 130), textcoords="offset points",
                        ha="center", va="center", fontsize=8.6, color="#37474F")
        ax.set_title(title, fontsize=12.5, pad=6)
        ax.text(0.5, -0.42, note, ha="center", va="center", fontsize=9.6,
                color=color)

    fig.suptitle("同一个单元、同一组坐标：节点顺序改变的是“形函数与物理节点的对应”",
                 fontsize=13, y=1.03)
    fig.tight_layout()
    save(fig, "fig09-node-ordering-detj.png")


# ==================================================== 图 10 等参映射与雅可比
def fig10():
    fig, ax = plt.subplots(figsize=(11.6, 5.0))
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-0.42, 6.90); ax.set_ylim(-1.30, 2.10)

    # 母单元
    ax.add_patch(Rectangle((0, 0), 1.5, 1.5, facecolor="#E8F1FB",
                           edgecolor="#37474F", lw=1.4, zorder=2))
    parent = {1: (0, 0), 2: (1.5, 0), 3: (1.5, 1.5), 4: (0, 1.5)}
    off = {1: (-7, -13), 2: (9, -13), 3: (9, 9), 4: (-7, 9)}
    align = {1: "right", 2: "left", 3: "left", 4: "right"}
    for k, (x, y) in parent.items():
        ax.plot([x], [y], "o", ms=8, color=BLUE, zorder=4)
        ax.annotate(str(k), (x, y), xytext=off[k], textcoords="offset points",
                    ha=align[k], fontsize=10.5, color=BLUE)
    ax.annotate("母单元（ξ, η）", (0.75, -0.30), ha="center", fontsize=10.5,
                color="#22303A")
    ax.annotate("−1 ≤ ξ, η ≤ 1", (0.75, -0.56), ha="center", fontsize=9.2,
                color=GRAY)
    arrow(ax, (1.52, 1.78), (1.52, 2.02), GRAY, lw=1.0, ms=8)
    ax.annotate("η", (1.52, 2.06), ha="center", fontsize=10, color=GRAY)
    arrow(ax, (1.55, 0.0), (1.80, 0.0), GRAY, lw=1.0, ms=8)
    ax.annotate("ξ", (1.84, 0.0), fontsize=10, color=GRAY, va="center")
    ax.add_patch(Rectangle((0.60, 0.60), 0.26, 0.26, facecolor="none",
                           edgecolor=BLUE, lw=1.2, ls="--", zorder=5))

    # 实际单元（任意四边形）
    phys = np.array([[3.50, 0.10], [5.55, -0.12], [6.25, 1.18], [3.72, 1.60]])
    ax.add_patch(Polygon(phys, closed=True, facecolor="#EAF4EC",
                         edgecolor="#37474F", lw=1.4, zorder=2))
    cen = phys.mean(axis=0)
    for k, (x, y) in enumerate(phys, start=1):
        ax.plot([x], [y], "o", ms=8, color=GREEN, zorder=4)
        d = np.array([x, y]) - cen
        d = d / np.linalg.norm(d)
        ax.annotate(str(k), (x, y), xytext=(9 * d[0], 9 * d[1]),
                    textcoords="offset points", ha="center", va="center",
                    fontsize=10.5, color=GREEN)
    ax.annotate("实际单元（x, y）", (4.85, -0.42), ha="center", fontsize=10.5,
                color="#22303A")
    ax.annotate("节点坐标 (xi, yi)", (4.85, -0.68), ha="center", fontsize=9.2,
                color=GRAY)
    ax.add_patch(Polygon([[4.45, 0.52], [4.63, 0.48], [4.68, 0.66],
                          [4.50, 0.70]], closed=True, facecolor="none",
                         edgecolor=GREEN, lw=1.2, ls="--", zorder=5))

    # 映射箭头
    arrow(ax, (1.72, 0.80), (3.30, 0.80), GRAY, lw=1.6, ms=12)
    arrow(ax, (3.30, 0.52), (1.72, 0.52), GRAY, lw=1.6, ms=12)
    ax.annotate("x = Σ Ni(ξ,η) xi\ny = Σ Ni(ξ,η) yi", (2.51, 1.52),
                ha="center", fontsize=10, color=GRAY, linespacing=1.6)
    ax.annotate("J = ∂(x,y)/∂(ξ,η)", (2.51, 0.06), ha="center",
                fontsize=11, color=RED)
    ax.annotate("dξ dη  →  det J · dξ dη", (2.51, -0.28), ha="center",
                fontsize=10, color=RED)

    fig.suptitle("等参映射：同一组形函数既插值坐标，也插值位移", fontsize=13, y=1.00)
    fig.tight_layout()
    save(fig, "fig10-isoparametric-jacobian.png")


# ==================================================== 图 11 组装：局部 → 全局
def fig11():
    fig = plt.figure(figsize=(12.4, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.22)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-0.55, 2.55); ax.set_ylim(-0.75, 1.75)

    nodes = {1: (0, 0), 2: (1, 0), 3: (2, 0), 4: (0, 1), 5: (1, 1), 6: (2, 1)}
    E1, E2 = [1, 2, 5, 4], [2, 3, 6, 5]
    for elem, color, ls in ((E1, BLUE, "-"), (E2, GREEN, "-")):
        ax.add_patch(Polygon([nodes[k] for k in elem], closed=True,
                             facecolor="none", edgecolor=color, lw=1.6,
                             ls=ls, zorder=2))
    ax.plot([1, 1], [0, 1], color="#37474F", lw=2.0, zorder=3)
    for k, (x, y) in nodes.items():
        ax.plot([x], [y], "o", ms=8, color="#37474F", zorder=4)
        ax.annotate("节点 %d\n(%d,%d)" % (k, 2 * k - 1, 2 * k), (x, y),
                    xytext=(0, 12 if y > 0 else -22),
                    textcoords="offset points", ha="center",
                    va="bottom" if y > 0 else "top",
                    fontsize=8.6, color="#37474F")
    ax.annotate("E1 → 全局 [1,2 | 3,4 | 9,10 | 7,8]", (0.5, 0.64),
                ha="center", va="center", fontsize=8.6, color=BLUE)
    ax.annotate("E2 → 全局 [3,4 | 5,6 | 11,12 | 9,10]", (1.5, 0.34),
                ha="center", va="center", fontsize=8.6, color=GREEN)
    ax.annotate("节点 2、5 被两个单元共享\n→ 公共自由度的贡献必须“累加”",
                (1.0, -0.42), ha="center", fontsize=9.6, color=RED)
    ax.set_title("两个 Q4 单元、六个节点、12 个自由度", fontsize=12)

    # 稀疏结构图
    ax = fig.add_subplot(gs[0, 1])
    n = 12
    S = np.zeros((n, n))
    for elem in (E1, E2):
        dofs = [(2 * k - 2, 2 * k - 1) for k in elem]
        flat = [d for pair in dofs for d in pair]
        for a in flat:
            for b in flat:
                S[a, b] = 1
    for i in range(n):
        for j in range(n):
            if S[i, j]:
                shared = (i in (2, 3, 8, 9)) and (j in (2, 3, 8, 9))
                ax.add_patch(Rectangle((j, n - 1 - i), 1, 1,
                                       facecolor=("#F2C6C6" if shared else "#D6E4F2"),
                                       edgecolor="white", lw=0.7))
    ax.set_xlim(0, n); ax.set_ylim(0, n)
    ax.set_aspect("equal")
    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_xticklabels(range(1, n + 1), fontsize=7.5)
    ax.set_yticklabels(range(n, 0, -1), fontsize=7.5)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("组装后的非零结构（红色 = 共享自由度的耦合块）", fontsize=11.5)
    ax.set_xlabel("全局自由度 J", fontsize=9.5)
    ax.set_ylabel("全局自由度 I", fontsize=9.5)
    fig.suptitle("组装不是“填表”，而是把每个单元的贡献按自由度映射累加到全局矩阵",
                 fontsize=13, y=0.99)
    fig.tight_layout()
    save(fig, "fig11-assembly-scatter.png")


# ==================================================== 图 12 三种稀疏存储
def fig12():
    fig = plt.figure(figsize=(13.6, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.25], wspace=0.20)

    # ------------------------------------------------ (a) Skyline
    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    n = 6
    starts = [0, 0, 1, 1, 1, 3]
    cnt = [i - s + 1 for i, s in enumerate(starts)]
    ptr = np.cumsum(cnt)          # 对角元在打包数组中的位置（1 起）
    k = 0
    for i in range(n):
        for j in range(n):
            if j > i:
                face, txt = "#F1F3F5", ""
            elif j < starts[i]:
                face, txt = "white", ""
            else:
                k += 1
                face, txt = "#D6E4F2", str(k)
            ax.add_patch(Rectangle((j, n - 1 - i), 1, 1, facecolor=face,
                                   edgecolor="#C9D2DA", lw=0.8))
            if txt:
                ax.text(j + 0.5, n - 1 - i + 0.5, txt, ha="center",
                        va="center", fontsize=8.6, color="#22303A")
    for j in range(n):
        ax.text(j + 0.5, n + 0.22, str(j + 1), ha="center", va="bottom",
                fontsize=8, color=GRAY)
    for i in range(n):
        ax.text(-0.18, n - 1 - i + 0.5, str(i + 1), ha="right", va="center",
                fontsize=8, color=GRAY)
    ax.set_xlim(-0.05, n + 0.05); ax.set_ylim(-1.45, n + 0.62)
    ax.set_aspect("equal")
    ax.set_title("① Skyline：每行从首个非零元\n一直存到对角元", fontsize=10.5)
    ax.text(n / 2, -0.55, "ptr = [%s]（对角元位置，ptr(0)=0）"
            % ", ".join(str(p) for p in ptr),
            ha="center", va="top", fontsize=9.0, color=BLUE)
    ax.text(n / 2, -1.00, "ptr(j) − ptr(j−1) = 第 j 行存储个数",
            ha="center", va="top", fontsize=8.4, color=GRAY)

    # ------------------------------------------------ (b) COO
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    pattern = [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (3, 1), (3, 3),
               (4, 2), (4, 4)]
    vals = list(range(1, 10))
    n = 5
    for i in range(n):
        for j in range(n):
            ax.add_patch(Rectangle((j, n - 1 - i), 1, 1, facecolor="white",
                                   edgecolor="#D5DCE3", lw=0.8))
    for (i, j), v in zip(pattern, vals):
        ax.add_patch(Rectangle((j, n - 1 - i), 1, 1, facecolor="#D6E4F2",
                               edgecolor="#37474F", lw=0.8))
        ax.text(j + 0.5, n - 1 - i + 0.5, str(v), ha="center", va="center",
                fontsize=9, color="#22303A")
    for j in range(n):
        ax.text(j + 0.5, n + 0.20, str(j + 1), ha="center", va="bottom",
                fontsize=8, color=GRAY)
    for i in range(n):
        ax.text(-0.18, n - 1 - i + 0.5, str(i + 1), ha="right", va="center",
                fontsize=8, color=GRAY)
    ax.set_xlim(-0.05, n + 0.05); ax.set_ylim(-2.95, n + 0.60)
    ax.set_aspect("equal")
    ax.set_title("② COO：每个非零元存\n“行、列、值”三个数", fontsize=10.5)
    ax.text(-0.05, -0.30, "ROW = [%s]" % ",".join(str(i + 1) for i, _ in pattern),
            fontsize=8.6, color=BLUE, family="monospace", ha="left", va="top")
    ax.text(-0.05, -0.90, "COL = [%s]" % ",".join(str(j + 1) for _, j in pattern),
            fontsize=8.6, color=GREEN, family="monospace", ha="left", va="top")
    ax.text(-0.05, -1.50, "VAL = [%s]" % ",".join(str(v) for v in vals),
            fontsize=8.6, color="#8D6E63", family="monospace", ha="left", va="top")
    ax.text(-0.05, -2.15, "好处：可乱序、可重复累加\n代价：row + col 各占一个整数",
            fontsize=8.2, color=GRAY, ha="left", va="top")

    # ------------------------------------------------ (c) CSR
    ax = fig.add_subplot(gs[0, 2])
    ax.axis("off")
    ax.set_xlim(-2.0, 9.6); ax.set_ylim(-1.55, 4.1)
    ax.set_title("③ CSR：按行压缩，VAL 与 COL 等长，\nPTR 记录每行起点", fontsize=10.5)
    val = list(range(1, 10))
    col = [1, 1, 2, 2, 3, 2, 4, 3, 5]
    ptr = [1, 2, 4, 6, 8]
    for r, (label, data, color) in enumerate(
            (("VAL", val, "#FDF3E3"), ("COL", col, "#EAF4EC"),
             ("PTR", ptr, "#EDE7F6"))):
        for c, txt in enumerate(data):
            ax.add_patch(Rectangle((c, 2.4 - r), 1, 1, facecolor=color,
                                   edgecolor="#C9D2DA", lw=0.7))
            ax.text(c + 0.5, 2.9 - r, txt, ha="center", va="center",
                    fontsize=8.8, color="#22303A")
        ax.text(-0.15, 2.9 - r, label, ha="right", va="center",
                fontsize=9.2, color="#37474F")
    ax.text(-0.15, 0.10,
            "查 K(3,2)：第 3 行 = COL(4…5) = [2, 3] → K(3,2) = VAL(4) = 4\n"
            "查 K(1,2)：第 1 行 = COL(1…1) = [1] 　→ 不存在 → 0",
            fontsize=8.6, color=RED, ha="left", va="top", linespacing=1.8)
    ax.text(-0.05, -0.75, "第 3 行：[PTR(3), PTR(4)−1] = [4, 5] → 2 个元素\n"
                          "第 5 行：[8, nnz]；完整 CSR 再补末位 PTR(6)=10（长 n+1）",
            fontsize=8.4, color=GRAY, ha="left", va="top", linespacing=1.8)

    fig.suptitle("三种稀疏存储：Skyline 存轮廓，COO 存三元组，CSR 存行压缩",
                 fontsize=13, y=1.00)
    fig.tight_layout()
    save(fig, "fig12-sparse-storage.png")


# ==================================================== 图 13 约束处理三法
def fig13():
    fig, axes = plt.subplots(1, 4, figsize=(13.0, 4.6))
    base = np.array([[2.0, -1.0, 0.0],
                     [-1.0, 2.0, -1.0],
                     [0.0, -1.0, 2.0]])
    titles = ["原方程  K u = F",
              "① 消元法\n删去受约束的行与列",
              "② 改行列法\n该行该列清零，对角置 1",
              "③ 大数（罚）法\n对角加 α，右端加 α·u0"]
    notes = ["约束自由度 i = 2 且 u₂ = 0.25",
             "未知量降为 2 个；其余自由度的 F 要减去 K·u0",
             "矩阵规模不变；F 同样要减 K·u0",
             "解近似满足；误差 ~ K/α，条件数变差"]
    for ax, title, note in zip(axes, titles, notes):
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_xlim(-0.9, 3.4); ax.set_ylim(-1.35, 3.9)
        ax.set_title(title, fontsize=10.5, pad=8)
        for i in range(3):
            for j in range(3):
                val = base[i, j]
                face = "white"
                txtc = "#22303A"
                alpha = 1.0
                if "消元" in title and i == 1:
                    face, alpha = "#E9EDF1", 0.6; txtc = "#AAB6C1"
                if "消元" in title and j == 1:
                    face, alpha = "#E9EDF1", 0.6; txtc = "#AAB6C1"
                if "改行列" in title and (i == 1 or j == 1):
                    face = "#E9EDF1"; txtc = "#AAB6C1"
                if "大数" in title and (i == 1 and j == 1):
                    face = "#D8EFD9"; txtc = GREEN
                ax.add_patch(Rectangle((j, 2 - i), 1, 1, facecolor=face,
                                       edgecolor="#C9D2DA", lw=0.8,
                                       alpha=alpha))
                show = val
                if "改行列" in title and (i == 1 or j == 1):
                    show = 1.0 if (i == 1 and j == 1) else 0.0
                if "大数" in title and (i == 1 and j == 1):
                    show = None
                if show is not None:
                    ax.text(j + 0.5, 2 - i + 0.5, "%.0f" % show,
                            ha="center", va="center", fontsize=9.6,
                            color=txtc)
                else:
                    ax.text(j + 0.5, 2 - i + 0.5, "α", ha="center",
                            va="center", fontsize=12, color=GREEN)
        ax.text(1.5, -0.55, note, ha="center", va="top", fontsize=8.8,
                color=RED if "误差" in note else GRAY, wrap=True)
        if "消元" in title:
            ax.plot([-0.12, 3.12], [1.0, 1.0], color=RED, lw=1.0,
                    ls=(0, (4, 3)))
            ax.plot([1.0, 1.0], [-0.12, 3.12], color=RED, lw=1.0,
                    ls=(0, (4, 3)))
    fig.suptitle("同一个约束，三种写法：改的是矩阵与右端项，不是物理", fontsize=13, y=1.03)
    fig.tight_layout()
    save(fig, "fig13-constraint-methods.png")


# ==================================================== 图 14 荷载向量
def fig14():
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 5.2),
                             gridspec_kw={"width_ratios": [1.05, 0.60, 1.0]})
    ax = axes[0]
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-0.55, 2.55); ax.set_ylim(-0.85, 1.75)
    ax.add_patch(Rectangle((0, 0), 2, 1, facecolor="#F4F7FA",
                           edgecolor="#37474F", lw=1.4, zorder=2))
    nodes = {1: (0, 0), 2: (2, 0), 3: (2, 1), 4: (0, 1)}
    for k, (x, y) in nodes.items():
        ax.plot([x], [y], "o", ms=8, color="#37474F", zorder=4)
        ax.annotate(str(k), (x, y), xytext=(-14, -14),
                    textcoords="offset points", fontsize=9.5, color="#37474F")
    # 体力
    for x in np.linspace(0.24, 1.76, 6):
        for y in (0.36, 0.64):
            arrow(ax, (x, y), (x, y - 0.20), "#8D6E63", lw=1.0, ms=7)
    ax.annotate("体力 b = 2000 N/m³（向下）\n总体 400 N → 每节点 100 N",
                xy=(0.55, 0.34), xytext=(-0.48, -0.62),
                fontsize=9, color="#6D4C41", ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color="#8D6E63", lw=1.0))
    # 右边界面力
    for y in np.linspace(0.12, 0.88, 5):
        arrow(ax, (2.0, y), (2.32, y), BLUE, lw=1.2, ms=8)
    ax.annotate("面力 t = 1000 N/m²\n合力 100 N → 节点 2、3 各 50 N",
                (2.42, 0.5), fontsize=9, color=BLUE, ha="left", va="center")
    # 节点集中力
    arrow(ax, (2.0, 1.0), (2.0, 1.30), RED, lw=2.0, ms=11)
    ax.annotate("+30 N", (2.0, 1.34), ha="center", fontsize=9.5, color=RED)
    ax.set_title("2 m × 1 m × 0.1 m 的 Q4：三种荷载同时作用", fontsize=12)

    # ---- 中：点荷载一致性
    ax = axes[1]
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-0.20, 1.16); ax.set_ylim(-0.02, 2.34)
    ax.add_patch(Rectangle((0.08, 1.55), 0.84, 0.70, facecolor="#F4F7FA",
                           edgecolor="#37474F", lw=1.2, zorder=2))
    ax.text(0.50, 2.12, "P = 30 N", ha="center", va="center",
            fontsize=10, color=RED, zorder=6)
    arrow(ax, (0.5, 2.02), (0.5, 1.78), RED, lw=2.2, ms=12)
    ax.plot([0.5], [1.86], "o", ms=4.5, color=RED, zorder=5)
    ax.text(-0.06, 1.90, "①", ha="center", va="center", fontsize=10.5,
            color="#22303A")

    ax.add_patch(Rectangle((0.08, 0.36), 0.84, 0.70, facecolor="#F4F7FA",
                           edgecolor="#37474F", lw=1.2, zorder=2))
    for x in (0.08, 0.92):
        ax.plot([x], [0.36], "o", ms=5, color="#37474F", zorder=5)
        ax.plot([x], [1.06], "o", ms=5, color="#37474F", zorder=5)
        arrow(ax, (x, 1.34), (x, 1.09), RED, lw=1.8, ms=10, zorder=6)
        arrow(ax, (x, 0.32), (x, 0.09), RED, lw=1.8, ms=10, zorder=6)
    ax.text(0.50, 0.71, "四节点\n各 7.5 N", ha="center", va="center",
            fontsize=9.5, color=RED, linespacing=1.5)
    ax.text(-0.06, 0.71, "②", ha="center", va="center", fontsize=10.5,
            color="#22303A")
    ax.text(0.50, 0.02, "$N_a(0,0)=1/4$，合力仍为 30 N",
            ha="center", va="bottom", fontsize=9.2, color="#4A5A66")
    ax.set_title("点荷载一致性", fontsize=12)

    ax = axes[2]
    labels = ["u1x", "u1y", "u2x", "u2y", "u3x", "u3y", "u4x", "u4y"]
    F = np.array([0, -100, 50, -100, 50, -70, 0, -100])
    colors = ["#8D6E63" if abs(v) == 100 and i % 2 == 1 else "#6D4C41"
              for i, v in enumerate(F)]
    for i, v in enumerate(F):
        c = "#8D6E63"
        if i in (2, 4):
            c = "#1565C0"
        elif i == 5:
            c = "#C62828"
        ax.bar(i, v, color=c, width=0.62, edgecolor="white")
        ax.annotate("%+d" % v, (i, v + (8 if v >= 0 else -8)),
                    ha="center", va="bottom" if v >= 0 else "top",
                    fontsize=9, color="#22303A")
    ax.axhline(0, color="#37474F", lw=1.0)
    ax.set_xticks(range(8)); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("节点荷载分量 / N", fontsize=10)
    ax.set_ylim(-125, 72)
    ax.grid(axis="y", color="#E4E9EE", lw=0.8)
    ax.set_axisbelow(True)
    ax.set_title("组装后的 F：ΣFx = +100 N，ΣFy = −370 N", fontsize=12)
    handles = [Rectangle((0, 0), 1, 1, facecolor=c) for c in
               ("#8D6E63", "#1565C0", "#C62828")]
    fig.legend(handles, ["体力（各 −100 N）", "边界面力（各 +50 N）",
                         "集中力（+30 N，与体力抵消）"],
               fontsize=8.8, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("三种荷载，三个公式：集中力找作用点，面力积边界，体力积单元",
                 fontsize=13, y=1.01)
    fig.tight_layout(rect=(0, 0.06, 1, 0.98))
    save(fig, "fig14-load-vector.png")


# ==================================================== 图 15 求解器代价量级
def fig15():
    N = np.logspace(2, 8, 200)
    dense = N ** 3
    skyline = N ** 2
    cg = 40 * N ** 1.5
    ref = dense[(np.abs(N - 1e4)).argmin()]
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    ax.loglog(N, dense / ref, color=RED, lw=2.2, label="稠密 LU / Cholesky：O(N³)")
    ax.loglog(N, skyline / ref, color="#EF9A3F", lw=2.2,
              label="带状 / Skyline 直接法：O(N·b²)，二维 b≈√N")
    ax.loglog(N, cg / ref, color=BLUE, lw=2.2,
              label="CG / PCG：O(N·nnz·k)，k 为迭代次数")
    ax.axvline(1e4, color=GRAY, lw=0.9, ls=(0, (4, 3)))
    ax.annotate("本库示例量级\n（约 10⁴ 自由度）", (1e4, 1.4e-8), xytext=(18, 4),
                textcoords="offset points", fontsize=9, color=GRAY,
                ha="left", va="bottom")
    ax.set_xlabel("自由度 N", fontsize=10.5)
    ax.set_ylabel("相对计算量（以 N = 10⁴ 归一）", fontsize=10.5)
    ax.grid(True, which="both", color="#E8EDF2", lw=0.7)
    ax.set_axisbelow(True)
    ax.legend(fontsize=9.4, loc="upper left", frameon=True, framealpha=0.95)
    ax.set_title("代价量级：为什么“几十万自由度”必须换算法", fontsize=12.5)
    ax.text(0.99, 0.04, "只比较量级，不含常数；实际选择还取决于条件数、内存与多工况复用",
            transform=ax.transAxes, ha="right", fontsize=8.6, color=GRAY)
    fig.tight_layout()
    save(fig, "fig15-solver-cost.png")


# ==================================================== 图 16 CG / PCG 收敛
def _lap_layers(n, ratio=8.0):
    """层状系数的二维拉普拉斯矩阵（Dirichlet 外边界），用于演示 CG/PCG。"""
    N = n * n
    A = np.zeros((N, N))

    def coef(i):
        return 1.0 if (i // 4) % 2 == 0 else ratio

    for i in range(n):
        for j in range(n):
            k = i * n + j
            s = 0.0
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < n and 0 <= jj < n:
                    w = 1.0 if di else coef(i)
                    A[k, ii * n + jj] -= w
                    s += w
            A[k, k] = s
    bnd = [i * n + j for i in range(n) for j in range(n)
           if i in (0, n - 1) or j in (0, n - 1)]
    for k in bnd:
        A[k, :] = 0.0
        A[:, k] = 0.0
        A[k, k] = 1.0
    return A, np.array(bnd)


def _cg(A, b, tol=1e-8, maxit=4000, M=None):
    x = np.zeros_like(b)
    r = b - A @ x
    z = r if M is None else M @ r
    p = z.copy()
    rz = r @ z
    b2 = np.linalg.norm(b)
    hist = [np.linalg.norm(r) / b2]
    for _ in range(maxit):
        Ap = A @ p
        alpha = rz / (p @ Ap)
        x += alpha * p
        r -= alpha * Ap
        hist.append(np.linalg.norm(r) / b2)
        if hist[-1] < tol:
            break
        z = r if M is None else M @ r
        rz_new = r @ z
        p = z + (rz_new / rz) * p
        rz = rz_new
    return x, np.array(hist)


def fig16():
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.9))
    n = 24
    A, bnd = _lap_layers(n)
    b = np.ones(n * n)
    b[bnd] = 0.0
    x0, h0 = _cg(A, b)
    x1, h1 = _cg(A, b, M=np.diag(1.0 / np.diag(A)))

    ax = axes[0]
    ax.semilogy(np.arange(len(h0)), h0, color=RED, lw=1.8, label="CG")
    ax.semilogy(np.arange(len(h1)), h1, color=BLUE, lw=1.8, label="PCG（Jacobi）")
    ax.axhline(1e-8, color=GRAY, lw=0.9, ls=(0, (4, 3)))
    ax.annotate("容差 1e−8", (0.98, 1e-8), xycoords=("axes fraction", "data"),
                xytext=(0, 7), textcoords="offset points",
                ha="right", fontsize=8.8, color=GRAY)
    ax.set_xlabel("迭代次数", fontsize=10)
    ax.set_ylabel("相对残差 ‖b − Ku‖ / ‖b‖", fontsize=10)
    ax.set_title("%d×%d 层状网格：CG %d 步、PCG %d 步" % (n, n, len(h0) - 1, len(h1) - 1),
                 fontsize=11)
    ax.grid(True, which="both", color="#E8EDF2", lw=0.7)
    ax.set_axisbelow(True)
    ax.legend(fontsize=9.5)

    ax = axes[1]
    ns = [16, 24, 32]
    it_cg, it_pcg, Ns = [], [], []
    for nn in ns:
        An, bn = _lap_layers(nn)
        bb = np.ones(nn * nn)
        bb[bn] = 0.0
        _, ha = _cg(An, bb)
        _, hb = _cg(An, bb, M=np.diag(1.0 / np.diag(An)))
        it_cg.append(len(ha) - 1)
        it_pcg.append(len(hb) - 1)
        Ns.append(nn * nn)
    ax.plot(Ns, it_cg, "o-", color=RED, lw=1.8, ms=6, label="CG")
    ax.plot(Ns, it_pcg, "s-", color=BLUE, lw=1.8, ms=5.5, label="PCG（Jacobi）")
    for xv, yv in zip(Ns, it_cg):
        ax.annotate(str(yv), (xv, yv), xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=8.6, color=RED)
    for xv, yv in zip(Ns, it_pcg):
        ax.annotate(str(yv), (xv, yv), xytext=(0, -14), textcoords="offset points",
                    ha="center", fontsize=8.6, color=BLUE)
    ax.set_xscale("log")
    ax.set_xticks(Ns)
    ax.set_xticklabels(["%d\n(%d×%d)" % (N, int(round(N ** 0.5)),
                                        int(round(N ** 0.5))) for N in Ns],
                       fontsize=8.8)
    ax.minorticks_off()
    ax.set_xlim(0.8 * min(Ns), 1.15 * max(Ns))
    ax.set_ylim(38, 146)
    ax.set_xlabel("自由度 N（括号内为网格边长）", fontsize=10)
    ax.set_ylabel("达到 1e−8 所需迭代次数", fontsize=10)
    ax.set_title("迭代次数随规模缓慢增长", fontsize=11)
    ax.grid(True, which="both", color="#E8EDF2", lw=0.7)
    ax.set_axisbelow(True)
    ax.legend(fontsize=9.5, loc="upper left")

    fig.suptitle("真实残差才是收敛判据：同一条曲线，预处理后少走约 1/4 的步数",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    save(fig, "fig16-cg-pcg-convergence.png")


# ============================================================ 图 17 雅可比：导数与映射
def _jac_model(x):
    """贯穿雅可比笔记的算例：y1 = x1^2 + x2, y2 = x1*x2"""
    return np.array([x[0] ** 2 + x[1], x[0] * x[1]])


def _jac_matrix(x):
    return np.array([[2.0 * x[0], 1.0], [x[1], x[0]]])


def fig17():
    """左：一元切线；右：小方格经 J 映射成平行四边形"""
    with plt.rc_context({"font.size": 13}):
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))

        # ---- 左：一元导数就是局部线性化
        ax = axes[0]
        f = lambda x: 0.45 * x ** 2 + 0.25 * x + 0.6
        df = lambda x: 0.9 * x + 0.25
        x0, dx = 1.0, 0.85
        xs = np.linspace(0.15, 2.35, 200)
        ax.plot(xs, f(xs), color=BLUE, lw=2.4, label="真实 $y=f(x)$")
        tx = np.array([x0 - 0.28, x0 + dx + 0.22])
        ax.plot(tx, f(x0) + df(x0) * (tx - x0), color=RED, lw=2.2, ls="--",
                label="切线（一阶近似）")
        ax.plot([x0], [f(x0)], "o", color=RED, ms=8, zorder=8)
        ax.annotate("$x_0$", (x0, f(x0)), xytext=(-6, -20), textcoords="offset points",
                    fontsize=14, color=GRAY, ha="center")
        ax.annotate("", (x0 + dx, f(x0)), (x0, f(x0)),
                    arrowprops=dict(arrowstyle="<|-|>", color=GRAY, lw=1.3,
                                    shrinkA=0, shrinkB=0))
        ax.annotate(r"$\Delta x$", (x0 + dx / 2, f(x0)), xytext=(0, -22),
                    textcoords="offset points", ha="center", fontsize=14, color=GRAY)
        ax.annotate("", (x0 + dx, f(x0 + dx)), (x0 + dx, f(x0)),
                    arrowprops=dict(arrowstyle="<|-|>", color=GREEN, lw=1.7,
                                    shrinkA=0, shrinkB=0))
        ax.annotate("真实\n$\\Delta y$", (x0 + dx, (f(x0) + f(x0 + dx)) / 2),
                    xytext=(10, 0), textcoords="offset points", fontsize=13,
                    color=GREEN, va="center")
        ax.set_xlim(0.1, 2.55)
        ax.set_ylim(0.0, 3.05)
        ax.set_xlabel("$x$", fontsize=14)
        ax.set_ylabel("$y$", fontsize=14)
        ax.set_title("一元：用一条切线代替一段曲线", fontsize=14)
        ax.legend(fontsize=11.5, loc="upper left", framealpha=0.95)
        ax.grid(True, color="#E8EDF2", lw=0.7)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=11)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

        # ---- 右：输入方格 -> 输出平行四边形（左右两套坐标框）
        ax = axes[1]
        x0v = np.array([2.0, 1.0])
        Jm = _jac_matrix(x0v)
        h = 1.0
        o_in = np.array([-4.3, 1.0])
        sq = np.array([[0, 0], [h, 0], [h, h], [0, h]]) + o_in
        poly = np.array([(0, 0), (h, 0), (h, h), (0, h)], float) @ Jm.T

        ax.add_patch(Polygon(sq, closed=True, facecolor=BLUE, alpha=0.18,
                             edgecolor=BLUE, lw=2.0, joinstyle="round"))
        ax.add_patch(Polygon(poly, closed=True, facecolor=RED, alpha=0.13,
                             edgecolor=RED, lw=2.0, joinstyle="round"))
        arrow(ax, o_in, o_in + np.array([1.9, 0]), GRAY, lw=1.1, ms=10)
        arrow(ax, o_in, o_in + np.array([0, 1.9]), GRAY, lw=1.1, ms=10)
        arrow(ax, (0, 0), (5.6, 0), GRAY, lw=1.1, ms=10)
        arrow(ax, (0, 0), (0, 4.4), GRAY, lw=1.1, ms=10)
        arrow(ax, (0, 0), poly[1], RED, lw=2.2, ms=13)
        arrow(ax, (0, 0), poly[3], RED, lw=2.2, ms=13)
        ax.annotate(r"$\Delta x_1$", (o_in[0] + h / 2, o_in[1] - 0.22), ha="center",
                    va="top", fontsize=13, color=BLUE)
        ax.annotate(r"$\Delta x_2$", (o_in[0] - 0.18, o_in[1] + h / 2), ha="right",
                    va="center", fontsize=13, color=BLUE)
        ax.annotate(r"$\mathbf{J}(\Delta x_1,0)^{\mathrm{T}}$",
                    poly[1] + np.array([0.00, -0.55]),
                    ha="center", va="top", fontsize=12, color=RED)
        ax.annotate(r"$\mathbf{J}(0,\Delta x_2)^{\mathrm{T}}$",
                    poly[3] + np.array([-0.10, 0.20]),
                    ha="left", va="bottom", fontsize=12, color=RED)
        arrow(ax, (-2.25, 1.55), (-0.75, 1.55), "#37474F", lw=2.4, ms=16)
        ax.annotate("$\\mathbf{J}$", (-1.50, 1.90), ha="center", fontsize=17,
                    color="#37474F")
        ax.annotate("输入 $(x_1,x_2)$", (-3.55, -0.55), ha="center",
                    fontsize=12, color=GRAY)
        ax.annotate("输出 $(y_1,y_2)$", (2.7, -0.55), ha="center",
                    fontsize=12, color=GRAY)
        ax.set_xlim(-5.1, 6.6)
        ax.set_ylim(-0.95, 4.6)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("二元：两列就是两个坐标方向的像", fontsize=14)
        for s in ("top", "right", "bottom", "left"):
            ax.spines[s].set_visible(False)

        fig.tight_layout()
        save(fig, "fig17-jacobian-tangent-and-map.png")


# ============================================================ 图 18 det J 的四种情形
def fig18():
    """同一个函数，四个工作点：放大 / 压缩 / 翻转 / 退化（同一坐标窗，直接可比）"""
    cases = [
        ((2.0, 1.0), "放大", "+7", GREEN),
        ((0.5, 0.0), "压缩", "+0.5", BLUE),
        ((0.0, 1.0), "翻转", "-1", "#E65100"),
        ((1.0, 2.0), "退化", "0", RED),
    ]
    with plt.rc_context({"font.size": 13}):
        fig, axes = plt.subplots(2, 2, figsize=(8.6, 8.2))
        h = 1.0
        sq = np.array([[0, 0], [h, 0], [h, h], [0, h]], float)
        lim = 5.3
        for k, (ax, (p, tag, dtxt, col)) in enumerate(zip(axes.ravel(), cases)):
            x0 = np.array(p, float)
            Jm = _jac_matrix(x0)
            img = (Jm @ sq.T).T
            ax.add_patch(Polygon(sq, closed=True, facecolor=BLUE, alpha=0.16,
                                 edgecolor=BLUE, lw=1.6, ls=(0, (4, 3))))
            ax.add_patch(Polygon(img, closed=True, facecolor=col, alpha=0.22,
                                 edgecolor=col, lw=2.1, joinstyle="round"))
            arrow(ax, (0, 0), tuple(img[1]), col, lw=1.9, ms=11)
            arrow(ax, (0, 0), tuple(img[3]), col, lw=1.9, ms=11)
            kmax = int(np.argmax(np.hypot(*(img - sq).T)))
            a0, a1 = sq[kmax], img[kmax]
            ax.plot([a0[0]], [a0[1]], "o", ms=6, color=BLUE, zorder=8)
            ax.plot([a1[0]], [a1[1]], "o", ms=7, color=col, zorder=9)
            ax.annotate("$A$", a0, xytext=(4, 5), textcoords="offset points",
                        fontsize=12, color=BLUE)
            ax.annotate("$A'$", a1, xytext=(-8, 6), textcoords="offset points",
                        fontsize=12.5, color=col, ha="right", va="bottom")
            if np.hypot(*(a1 - a0)) > 0.2:
                ax.annotate("", a1, a0,
                            arrowprops=dict(arrowstyle="-|>", color=col, lw=1.2,
                                            ls=(0, (3, 2)), shrinkA=3, shrinkB=3,
                                            alpha=0.8), zorder=7)
            ax.axhline(0, color=GRAY, lw=0.9)
            ax.axvline(0, color=GRAY, lw=0.9)
            ax.annotate("$\\det\\mathbf{J}=%s$" % dtxt, (0.5, 0.975),
                        xycoords="axes fraction", ha="center", va="top",
                        fontsize=14, color=col,
                        bbox=dict(boxstyle="round,pad=0.3", fc="white",
                                  ec=col, lw=1.2))
            ax.annotate(tag, (0.5, 0.875), xycoords="axes fraction",
                        ha="center", va="top", fontsize=14, color="#263238")
            ax.annotate("工作点 $(%g,\\ %g)$" % (p[0], p[1]), (0.5, 0.80),
                        xycoords="axes fraction", ha="center", va="top",
                        fontsize=11.5, color=GRAY)
            ax.set_xlim(-0.9, lim)
            ax.set_ylim(-0.9, lim)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ("top", "right", "bottom", "left"):
                ax.spines[s].set_visible(False)

        fig.suptitle("同一个函数： $\\det\\mathbf{J}=2x_1^2-x_2$ 决定局部面积怎么变",
                     fontsize=15, y=0.985)
        fig.text(0.5, 0.018,
                 "蓝色虚线为单位方格（输入），实色为其像（输出）；\n"
                 "同一坐标窗，大小可直接比较； $A\\to A'$ 看定向是否翻转。",
                 ha="center", fontsize=11.5, color=GRAY, linespacing=1.5)
        fig.tight_layout(rect=(0, 0.055, 1, 0.955))
        save(fig, "fig18-jacobian-det-cases.png")


# ============================================================ 图 19 二连杆奇异位形
def fig19():
    l1 = l2 = 1.0

    def pose(t1, t2):
        p0 = np.array([0.0, 0.0])
        p1 = p0 + l1 * np.array([np.cos(t1), np.sin(t1)])
        p2 = p1 + l2 * np.array([np.cos(t1 + t2), np.sin(t1 + t2)])
        return p0, p1, p2

    def jac(t1, t2):
        s1, c1 = np.sin(t1), np.cos(t1)
        s12, c12 = np.sin(t1 + t2), np.cos(t1 + t2)
        return np.array([[-l1 * s1 - l2 * s12, -l2 * s12],
                         [l1 * c1 + l2 * c12, l2 * c12]])

    with plt.rc_context({"font.size": 13}):
        fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.2))
        poses = [(np.radians(18), np.radians(62), "一般姿态：两连杆不共线"),
                 (np.radians(28), 0.0, "奇异姿态：两连杆伸直")]
        for ax, (t1, t2, title) in zip(axes, poses):
            p0, p1, p2 = pose(t1, t2)
            for a, b in ((p0, p1), (p1, p2)):
                ax.plot([a[0], b[0]], [a[1], b[1]], color="#455A64", lw=10,
                        solid_capstyle="round", zorder=3)
                ax.plot([a[0], b[0]], [a[1], b[1]], color="#90A4AE", lw=3.4,
                        solid_capstyle="round", zorder=4)
            for p in (p0, p1, p2):
                ax.plot([p[0]], [p[1]], "o", ms=10, color=BLUE, mec="white",
                        mew=1.5, zorder=6)
            Jm = jac(t1, t2)
            detJ = np.linalg.det(Jm)
            # 关节速度 |q_i| <= 1 时，末端可达速度集合 = J 作用在单位方格上
            sqv = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
            vel = (Jm @ (sqv.T * 0.34)).T
            col = RED if abs(detJ) < 1e-6 else GREEN
            if abs(detJ) < 1e-6:
                # 秩 1：可达速度集合退化成一条线段，直接画线而不画“零面积四边形”
                ax.plot([p2[0], (p2 + vel[2])[0]], [p2[1], (p2 + vel[2])[1]],
                        color=col, lw=5.2, alpha=0.4, solid_capstyle="round",
                        zorder=5)
                arrow(ax, p2, p2 + vel[2], col, lw=1.9, ms=11)
            else:
                ax.add_patch(Polygon(p2 + vel, closed=True, facecolor=col,
                                     alpha=0.24, edgecolor=col, lw=1.7,
                                     joinstyle="round", zorder=5))
                arrow(ax, p2, p2 + vel[1], col, lw=1.7, ms=10)
                arrow(ax, p2, p2 + vel[3], col, lw=1.7, ms=10)
            ax.annotate("关节 1", p0, xytext=(8, -20), textcoords="offset points",
                        fontsize=11.5, color=GRAY)
            ax.annotate("关节 2", p1, xytext=(-10, 12), textcoords="offset points",
                        fontsize=11.5, color=GRAY, ha="right")
            ax.annotate("末端", p2, xytext=(12, -18), textcoords="offset points",
                        fontsize=11.5, color=GRAY)
            if abs(detJ) < 1e-6:
                seg_mid = p2 + vel[2]
                note = ("$\\dot q_1,\\dot q_2$ 任取：末端速度\n"
                        "被压成一条线段（秩 1），\n"
                        "只能垂直于连杆运动")
            else:
                seg_mid = p2 + vel[2]
                note = ("$|\\dot q_i|\\leq 1$ 时末端速度\n"
                        "可达集是一个平行四边形\n"
                        "（秩 2，面积 $=|\\det\\mathbf{J}|$）")
            ax.annotate(note, xy=seg_mid, xycoords="data",
                        xytext=(0.025, 0.975), textcoords="axes fraction",
                        ha="left", va="top", fontsize=11.8, color=col,
                        linespacing=1.55,
                        arrowprops=dict(arrowstyle="-", lw=1.1, color=col,
                                        shrinkA=4, shrinkB=9,
                                        connectionstyle="arc3,rad=-0.14"))
            ax.annotate("$\\det\\mathbf{J}=l_1l_2\\sin\\theta_2=%s$"
                        % ("0" if abs(detJ) < 1e-6 else "%.2f" % detJ),
                        (0.5, 0.035), xycoords="axes fraction", ha="center",
                        fontsize=13.5, color=col,
                        bbox=dict(boxstyle="round,pad=0.32", fc="white",
                                  ec=col, lw=1.1))
            ax.set_title(title, fontsize=13.5)
            ax.set_xlim(-1.75, 3.05)
            ax.set_ylim(-1.05, 2.25)
            ax.set_aspect("equal")
            ax.grid(True, color="#EEF2F5", lw=0.7)
            ax.set_axisbelow(True)
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)

        fig.suptitle("奇异位形： $\\det\\mathbf{J}=0$ 时，关节转得再快也伸不出连杆",
                     fontsize=14.5, y=0.995)
        fig.tight_layout()
        save(fig, "fig19-jacobian-singular-arm.png")


# ================================== 图 20 Ritz 与 Galerkin：两条路，一个方程组
def fig20():
    fig, ax = plt.subplots(figsize=(12.8, 6.0))
    ax.set_xlim(0, 12.8)
    ax.set_ylim(0, 6.0)
    ax.axis("off")

    def box(x, y, w, h, title, body, fc):
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.10,rounding_size=0.14",
            facecolor=fc, edgecolor="#37474F", lw=1.0, zorder=2))
        ax.text(x + w / 2, y + h - 0.33, title, ha="center", va="center",
                fontsize=12.4, color="#22303A", zorder=3)
        ax.text(x + w / 2, y + h / 2 - 0.34, body, ha="center", va="center",
                fontsize=10.0, color="#4A5A66", linespacing=1.7, zorder=3)

    w, h = 2.85, 1.62
    xs = [0.75, 4.25, 7.75]
    y_top, y_bot = 3.72, 1.32

    rows = [
        (y_top, "#E8F1FB", [
            ("Ritz：先有能量", "总势能 $J(v)=\\frac{1}{2}A(v,v)-F(v)$\n求 $J$ 在 $H$ 上的极小"),
            ("限制到子空间", "$J(\\mathbf{c})=\\frac{1}{2}\\mathbf{c}^{\\top}K\\mathbf{c}-\\mathbf{f}^{\\top}\\mathbf{c}$\n$\\mathbf{c}=(c_1,\\dots,c_N)$"),
            ("极值条件", "$\\partial J/\\partial c_i=0$\n二次型正定 → 极小点唯一"),
        ]),
        (y_bot, "#EAF4EC", [
            ("Galerkin：只用弱形式", "$A(u,v)=F(v)$，$\\forall v\\in H$\n不要求 $A$ 对称"),
            ("测试函数取在子空间", "$A(u_N,v)=F(v)$，$\\forall v\\in S_N$\n$u_N=\\sum_j c_j\\varphi_j$"),
            ("$c_i$ 任意 → 括号为零", "$\\sum_i c_i[\\sum_j A(\\varphi_j,\\varphi_i)c_j-F(\\varphi_i)]=0$"),
        ]),
    ]
    for y, color, items in rows:
        for x, (t, b) in zip(xs, items):
            box(x, y, w, h, t, b, color)
        for x in xs[:-1]:
            arrow(ax, (x + w + 0.10, y + h / 2), (x + w + 0.42, y + h / 2),
                  GRAY, lw=1.6, ms=11)

    ax.add_patch(FancyBboxPatch(
        (11.02, 2.62), 1.52, 1.42, boxstyle="round,pad=0.10,rounding_size=0.16",
        facecolor="#FDF3E3", edgecolor="#B7202A", lw=1.4, zorder=4))
    ax.text(11.78, 3.52, "同一个", ha="center", va="center", fontsize=10.5,
            color="#8A2A2A", zorder=5)
    ax.text(11.78, 3.16, "$K\\mathbf{c}=\\mathbf{f}$", ha="center", va="center",
            fontsize=14.5, color="#B7202A", zorder=5)
    ax.text(11.78, 2.85, "对称时完全一致", ha="center", va="center",
            fontsize=8.8, color="#8A2A2A", zorder=5)
    arrow(ax, (10.72, 4.53), (11.10, 3.86), GRAY, lw=1.5, ms=11)
    arrow(ax, (10.72, 2.13), (11.10, 2.80), GRAY, lw=1.5, ms=11)

    ax.text(6.4, 5.62, "两条路径：能量极小化 vs 加权余量投影",
            ha="center", fontsize=15)
    ax.text(6.4, 0.86,
            "上排只有在 $A$ 对称正定时才存在（否则没有能量泛函可极小化）；"
            "下排只要求 $A$ 有界、强制，非对称问题照样能导出方程组。",
            ha="center", fontsize=10.8, color=GRAY)
    ax.text(6.4, 0.44,
            "所以：Ritz 是 Galerkin 在“对称正定”情形下的特例，而不是另一套方法。",
            ha="center", fontsize=10.8, color="#B7202A")
    save(fig, "fig20-ritz-galerkin-two-roads.png")


# ================================== 图 21 能量正交投影：误差垂直于子空间
def fig21():
    fig, ax = plt.subplots(figsize=(11.6, 6.2))
    ax.set_xlim(0, 11.6)
    ax.set_ylim(0, 6.2)
    ax.axis("off")

    y0 = 2.30                      # 子空间（画成一条带）
    xO, xN, xv = 1.30, 5.20, 8.30
    u_top = 4.55

    ax.add_patch(Rectangle((0.55, y0 - 0.30), 10.2, 0.60, facecolor=BLUE,
                           alpha=0.10, edgecolor=BLUE, lw=1.1, zorder=1))
    ax.text(10.55, y0 + 0.42, "子空间 $S_N$", color=BLUE, fontsize=14,
            ha="right", va="bottom")

    # 三个向量：u、u_N、以及任取的 v
    arrow(ax, (xO, y0), (xN, u_top), RED, lw=3.0, ms=15)
    arrow(ax, (xO, y0), (xN, y0), GREEN, lw=3.0, ms=15)
    arrow(ax, (xO, y0), (xv, y0), "#37474F", lw=1.6, ms=12)
    arrow(ax, (xN, y0), (xN, u_top), "#8E24AA", lw=2.4, ms=14)

    # 勾股三角形 u - v - u_N
    ax.plot([xN, xv], [u_top, y0], color=GRAY, lw=1.3, ls=(0, (4, 3)), zorder=3)

    # 直角标记（在 u_N 处，竖直边与子空间法向）
    s = 0.26
    ax.add_patch(Rectangle((xN, y0), s, s, facecolor="white",
                           edgecolor="#37474F", lw=1.1, zorder=5))

    for (xx, yy), col in [((xO, y0), "#37474F"), ((xN, y0), GREEN),
                          ((xN, u_top), RED), ((xv, y0), "#37474F")]:
        ax.plot([xx], [yy], "o", ms=8.5, color=col, mec="white", mew=1.4,
                zorder=6)

    ax.text(xO - 0.18, y0 - 0.30, "$O$", fontsize=14, ha="right", va="top")
    ax.text(xN - 0.20, y0 - 0.30, "$u_N$", color=GREEN, fontsize=16,
            ha="right", va="top")
    ax.text(xN - 0.22, u_top + 0.05, "$u$", color=RED, fontsize=16,
            ha="right", va="bottom")
    ax.text(xN + 0.20, 0.5 * (y0 + u_top), "$u-u_N$", color="#8E24AA",
            fontsize=15, ha="left", va="center")
    ax.text(xv + 0.10, y0 - 0.32, "$v$（任取）", color="#37474F", fontsize=14,
            ha="left", va="top")

    ax.text(0.55, 5.90, r"$A(u-u_N,\,v)=0$ 对所有 $v\in S_N$",
            fontsize=17, color="#22303A")
    ax.text(0.55, 5.42,
            "不是“误差很小”，而是误差与整个子空间精确正交",
            fontsize=12.5, color=GRAY)

    ax.text(0.55, 1.22,
            r"于是 $\Vert u-v\Vert_E^2=\Vert u-u_N\Vert_E^2"
            r"+\Vert u_N-v\Vert_E^2\ \geq\ \Vert u-u_N\Vert_E^2$。",
            fontsize=13.2, color="#22303A")
    ax.text(0.55, 0.74,
            "在 $S_N$ 能表示的所有函数里， $u_N$ 离真解最近——这就是“最佳逼近”，"
            "勾股定理的直接推论。",
            fontsize=12.5, color=GRAY)
    save(fig, "fig21-energy-orthogonal-projection.png")


# ================================== 图 22 收敛：能量单调上升，误差按阶下降
def fig22():
    ns = np.array([3, 7, 15, 31, 63, 127])
    hs = 1.0 / (ns + 1)
    A = np.array([4.686291501, 4.871709919, 4.918968217,
                  4.930839888, 4.933811384, 4.934554481])
    eA = np.array([4.985e-01, 2.512e-01, 1.258e-01,
                   6.295e-02, 3.148e-02, 1.574e-02])
    eL2 = np.array([3.928e-02, 9.921e-03, 2.487e-03,
                    6.220e-04, 1.555e-04, 3.888e-05])
    exact = np.pi ** 2 / 2

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.0))

    ax = axes[0]
    ax.axhline(exact, color=GRAY, lw=1.3, ls="--")
    ax.plot(ns, A, "o-", color=GREEN, lw=2.2, ms=7,
            label="$A(u_N,u_N)$（离散应变能）")
    ax.text(3.15, exact + 0.008, "$A(u,u)=\\pi^2/2=%.6f$" % exact,
            fontsize=11.2, color=GRAY, va="bottom")
    ax.text(0.97, 0.13,
            "曲线与虚线的差距\n$=A(u,u)-A(u_N,u_N)=\\Vert u-u_N\\Vert_E^2\\geq0$",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
            color="#B7202A", linespacing=1.6)
    ax.set_xscale("log")
    ax.set_xlabel("内部节点数 $N$（网格越来越细 →）", fontsize=12)
    ax.set_ylabel("应变能", fontsize=12)
    ax.set_title("(a) 离散能量从下方单调逼近真值", fontsize=13)
    ax.set_ylim(4.66, 4.96)
    ax.legend(loc="lower right", fontsize=11)
    ax.grid(True, color="#EEF2F5", lw=0.8)
    ax.set_axisbelow(True)

    ax = axes[1]
    ax.loglog(hs, eA, "o-", color=BLUE, lw=2.2, ms=7,
              label="能量范数误差 $\\Vert u-u_N\\Vert_E$")
    ax.loglog(hs, eL2, "s-", color=RED, lw=2.2, ms=6.5,
              label="$L^2$ 误差 $\\Vert u-u_N\\Vert_{L^2}$")
    hh = np.array([hs[-1], 0.30])
    ax.loglog(hh, eA[-1] * (hh / hs[-1]) ** 1.0, "--", color=BLUE, lw=1.2,
              alpha=0.8)
    ax.loglog(hh, eL2[-1] * (hh / hs[-1]) ** 2.0, "--", color=RED, lw=1.2,
              alpha=0.8)
    ax.text(0.058, 0.135, "斜率 1", color=BLUE, fontsize=11.5)
    ax.text(0.052, 0.0033, "斜率 2", color=RED, fontsize=11.5)
    ax.set_xlabel("单元尺寸 $h$（→ 越来越小）", fontsize=12)
    ax.set_ylabel("误差", fontsize=12)
    ax.set_title("(b) 双对数斜率 = 收敛阶", fontsize=13)
    ax.legend(loc="lower left", fontsize=11)
    ax.grid(True, which="both", color="#EEF2F5", lw=0.8)
    ax.set_axisbelow(True)
    ax.invert_xaxis()

    fig.suptitle("一维模型问题 $-u''=\\pi^2\\sin\\pi x$：分片线性 Galerkin 的收敛",
                 fontsize=15, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "fig22-galerkin-convergence.png")


# ================================== 图 23 全域基 vs 分片基：稠密与带状
def fig23():
    n = 40
    i = np.arange(1, n + 1, dtype=float)
    Kp = np.outer(i, i) / (i[:, None] + i[None, :] - 1.0)      # 多项式基
    h = 1.0 / (n + 1)
    Kh = (2.0 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)) / h  # hat 基

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.3))

    for ax, M, title, cmap in [
            (axes[0], Kp, "(a) 全域多项式基：稠密、病态", "Reds"),
            (axes[1], Kh, "(b) 分片 hat 基：三对角、稀疏", "Blues")]:
        Mm = np.abs(M)
        Mm = Mm / Mm.max()
        ax.imshow(Mm ** 0.35, cmap=cmap, interpolation="nearest",
                  vmin=0, vmax=1)
        nnz = int(np.count_nonzero(np.abs(M) > 1e-12 *
                                   (np.abs(M).max() if np.abs(M).max() else 1)))
        dens = 100.0 * nnz / n ** 2
        cond = np.linalg.cond(M)
        ax.set_title(title, fontsize=13)
        ax.set_xticks([0, n - 1])
        ax.set_yticks([0, n - 1])
        ax.set_xticklabels([1, n], fontsize=9)
        ax.set_yticklabels([1, n], fontsize=9)
        ax.set_xlabel("未知量编号 $j$", fontsize=11.5)
        ax.set_ylabel("未知量编号 $i$", fontsize=11.5)
        mant, expo = ("%.1f" % (cond / 10 ** np.floor(np.log10(cond))),
                      int(np.floor(np.log10(cond))))
        ax.text(0.5, 0.94,
                "非零元 %d / %d（%.1f%%）\n条件数 $\\kappa\\approx%s\\times10^{%d}$"
                % (nnz, n * n, dens, mant, expo),
                transform=ax.transAxes, ha="center", va="top", fontsize=10.8,
                color="#22303A",
                bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="#C9D2DA",
                          alpha=0.94))

    fig.suptitle("基函数选在哪里，决定刚度矩阵长什么样（$N=%d$）" % n,
                 fontsize=15, y=0.99)
    fig.text(0.5, 0.015,
             "左：基函数铺满全场的代价 —— 任意两基函数重叠，矩阵全满；"
             "右：基函数只在自己的节点邻域非零，重叠只能是相邻节点。",
             ha="center", fontsize=10.8, color=GRAY)
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    save(fig, "fig23-dense-vs-banded-stiffness.png")


# ============================================ 图 24 可压缩控制体收支
def fig24():
    fig, ax = plt.subplots(figsize=(11.2, 5.4))
    ax.set_xlim(0, 11.2)
    ax.set_ylim(0, 5.4)
    ax.axis("off")

    cv = FancyBboxPatch((3.7, 1.35), 4.4, 2.9,
                        boxstyle="round,pad=0.06", fc="#E3F2FD",
                        ec=BLUE, lw=2.0)
    ax.add_patch(cv)
    ax.text(5.9, 3.55, "固定控制体 $V_i$", fontsize=15, ha="center",
            color="#0D2A4A")
    ax.text(5.9, 2.9, "保存的是平均值", fontsize=12.5, ha="center",
            color=GRAY)
    ax.text(5.9, 2.25,
            r"$\bar{\mathbf{U}}_i=(\rho,\ \rho u,\ \mathcal{E})^{\mathsf{T}}$",
            fontsize=15, ha="center")

    rows = [(4.35, r"质量：$\rho u$", RED),
            (3.05, r"动量：$\rho u^2+p$", "#6A1B9A"),
            (1.8, r"总能量：$(\mathcal{E}+p)u$", GREEN)]
    for y, lab, col in rows:
        arrow(ax, (1.55, y), (3.62, y), col, lw=2.6)
        arrow(ax, (8.1, y), (10.2, y), col, lw=2.6)
        ax.text(1.42, y, lab, fontsize=12.5, ha="right", va="center",
                color=col)
        ax.text(9.15, y + 0.22, "流出", fontsize=10.5, ha="center",
                color=col)

    ax.text(0.85, 5.0, "单位时间、单位面积的交换 $\\mathbf{F}_n$",
            fontsize=12.5, ha="center", color="#22303A")
    ax.text(10.4, 5.0, "五分量：$1+3+1$", fontsize=12, ha="center",
            color="#22303A")

    ax.text(5.9, 0.66,
            r"$\frac{\mathrm{d}\bar{\mathbf{U}}_i}{\mathrm{d}t}"
            r"=-\frac{1}{V_i}\sum_f\widehat{\mathbf{F}}_{n,f}S_f$",
            fontsize=18, ha="center")
    ax.text(5.9, 0.12, "面收支汇总 → 守恒量变化率", fontsize=11.5,
            ha="center", color=GRAY)

    fig.suptitle("可压缩有限体积：单元里存多少，面上就交换多少",
                 fontsize=15.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, "fig24-compressible-control-volume.png")


# ============================================ 图 25 N-S 模型封闭关系
def fig25():
    fig, ax = plt.subplots(figsize=(11.4, 5.8))
    ax.set_xlim(0, 11.4)
    ax.set_ylim(0, 5.8)
    ax.axis("off")

    def box(x, y, w, h, title, body, fc, ec):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0.08", fc=fc,
                                    ec=ec, lw=1.6))
        ax.text(x + w / 2, y + h - 0.28, title, fontsize=12.5,
                ha="center", va="top", color=ec, weight="bold")
        ax.text(x + w / 2, y + h / 2 - 0.15, body, fontsize=11.2,
                ha="center", va="center", color="#22303A")

    box(0.35, 3.4, 2.7, 1.85, "质量守恒",
        r"$\partial_t\rho+\nabla\cdot(\rho\mathbf{u})=0$",
        "#E3F2FD", BLUE)
    box(0.35, 1.85, 2.7, 1.4, "动量守恒（3 条）",
        r"$\partial_t(\rho\mathbf{u})+\nabla\cdot(\rho\mathbf{u}\otimes\mathbf{u})$"
        "\n" r"$=-\nabla p+\nabla\cdot\boldsymbol{\tau}$",
        "#E3F2FD", BLUE)
    box(0.35, 0.3, 2.7, 1.4, "总能量守恒",
        r"$\partial_t\mathcal{E}+\nabla\cdot[(\mathcal{E}+p)\mathbf{u}]$"
        "\n" r"$=\nabla\cdot(\boldsymbol{\tau}\cdot\mathbf{u})-\nabla\cdot\mathbf{q}$",
        "#E3F2FD", BLUE)

    box(4.15, 3.55, 3.0, 1.7, "热力学封闭",
        r"$p=\rho RT$" "\n" r"$e=c_vT$",
        "#FFF3E0", "#EF6C00")
    box(4.15, 0.4, 3.0, 1.55, "应力与热流",
        r"$\boldsymbol{\tau}=\mu[\nabla\mathbf{u}+(\nabla\mathbf{u})^{\mathsf{T}}"
        r"-\frac{2}{3}(\nabla\cdot\mathbf{u})\mathbf{I}]$"
        "\n" r"$\mathbf{q}=-\kappa\nabla T$",
        "#E8F5E9", GREEN)

    ax.add_patch(FancyBboxPatch((8.0, 2.0), 3.0, 1.7,
                                boxstyle="round,pad=0.08", fc="#F3E5F5",
                                ec="#6A1B9A", lw=2.0))
    ax.text(9.5, 3.25, "封闭的黏性 N–S", fontsize=14, ha="center",
            color="#6A1B9A", weight="bold")
    ax.text(9.5, 2.6, "7 个标量场 = 5 条守恒\n+ 2 条热力学关系", fontsize=11,
            ha="center", color="#22303A")

    arrow(ax, (3.12, 3.9), (4.08, 4.25), "#EF6C00", lw=2)
    arrow(ax, (3.12, 1.35), (4.08, 1.2), GREEN, lw=2)
    arrow(ax, (7.22, 4.2), (7.95, 3.15), "#EF6C00", lw=2)
    arrow(ax, (7.22, 1.2), (7.95, 2.45), GREEN, lw=2)

    ax.text(9.5, 0.8, "黏度 $\\mu$、导热系数 $\\kappa$\n由物性给定，不再另立方程",
            fontsize=10.8, ha="center", color=GRAY)

    fig.suptitle("封闭 = 每一个出现的量都有关系可以确定",
                 fontsize=15.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, "fig25-ns-closure-map.png")


# ============================================ 图 26 边界特征方向
def fig26():
    cases = [
        ("超声速入口：3 条进入", -2.0, 1.0, "inlet"),
        ("亚声速入口：2 条进入", -0.5, 1.0, "inlet"),
        ("亚声速出口：1 条进入", 0.5, 1.0, "outlet"),
        ("超声速出口：0 条进入", 2.0, 1.0, "outlet"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 6.2))
    for ax, (title, un, c, kind) in zip(axes.ravel(), cases):
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-0.34, 1.15)
        ax.axis("off")
        dom = (0, 1.15) if kind == "inlet" else (-1.15, 0)
        ax.add_patch(Rectangle((dom[0], -0.06), dom[1] - dom[0],
                               1.18, fc="#E3F2FD", ec="none", zorder=0))
        ax.plot([0, 0], [-0.06, 1.13], color="#22303A", lw=1.4, ls="--")
        ax.text(0, 1.08, "边界", fontsize=10.5, ha="center",
                color="#22303A")
        ax.text(dom[0] + (dom[1] - dom[0]) / 2, 1.08, "计算域",
                fontsize=10.5, ha="center", color=BLUE)
        speeds = [un - c, un, un + c]
        colors = [RED, "#6A1B9A", GREEN]
        labels = [r"$u_n-c$", r"$u_n$", r"$u_n+c$"]
        # 外法向在左侧入口指向 -x、右侧出口指向 +x，
        # 物理传播方向 = 特征速度 s 乘以外法向
        orient = -1.0 if kind == "inlet" else 1.0
        for k, (s, col, lab) in enumerate(zip(speeds, colors, labels)):
            y = 0.78 - 0.26 * k
            dx = orient * np.sign(s) if s else 0.0
            ax.annotate("", xy=(0.38 * dx, y), xytext=(0, y),
                        arrowprops=dict(arrowstyle="-|>", color=col,
                                        lw=2.2, mutation_scale=13))
            ax.text(0.45 * dx if dx else 0.05, y,
                    "%s = %.1f" % (lab, s), fontsize=10.2, color=col,
                    va="center", ha="left" if dx >= 0 else "right")
        # 外法向为正，s<0 即朝区域内部传播，入口与出口同判据
        nin = sum(1 for s in speeds if s < 0)
        ax.set_title(title, fontsize=12.5, color="#22303A")
        note = ("进入 $%d$ 条 → 需给定 $%d$ 条外部关系" % (nin, nin)
                if nin else "进入 $0$ 条 → 无需外部条件，由内点外推")
        ax.text(0, -0.21, note, fontsize=10.5, ha="center",
                va="center", color=RED)

    fig.suptitle("边界给几个条件，由信息传播方向决定（一维无黏 Euler）",
                 fontsize=15, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "fig26-characteristic-boundaries.png")


# ============================================ 图 27 亚声速出口构造
def fig27():
    fig, ax = plt.subplots(figsize=(12.4, 5.6))
    ax.set_xlim(0, 12.4)
    ax.set_ylim(0, 5.6)
    ax.axis("off")

    def box(x, y, w, h, text, fc, ec, fs=11.2):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0.07", fc=fc,
                                    ec=ec, lw=1.5))
        ax.text(x + w / 2, y + h / 2, text, fontsize=fs,
                ha="center", va="center", color="#22303A")

    box(0.3, 3.75, 2.5, 1.1,
        "内侧状态\n" r"$\rho_{\mathrm{in}},\ u_{\mathrm{in}},\ p_{\mathrm{in}}$",
        "#E3F2FD", BLUE)
    box(0.3, 1.55, 2.5, 1.0, "外部给定\n" r"$p_b$", "#FFF3E0", "#EF6C00")
    box(3.35, 3.75, 2.6, 1.1,
        "熵关系\n" r"$\rho_b=\rho_{\mathrm{in}}(p_b/p_{\mathrm{in}})^{1/\gamma}$",
        "#E8F5E9", GREEN)
    box(3.35, 1.55, 2.6, 1.0,
        "声速\n" r"$c_b=\sqrt{\gamma p_b/\rho_b}$",
        "#E8F5E9", GREEN)
    box(6.45, 3.75, 2.7, 1.1,
        "声学组合\n" r"$J_+=u+2c/(\gamma-1)$",
        "#F3E5F5", "#6A1B9A")
    box(6.45, 1.55, 2.7, 1.0,
        "边界速度\n" r"$u_b=u_{\mathrm{in}}+\frac{2(c_{\mathrm{in}}-c_b)}{\gamma-1}$",
        "#F3E5F5", "#6A1B9A")
    box(9.65, 2.55, 2.45, 1.35,
        "守恒变量\n" r"$\mathbf{U}_b=(\rho_b,\rho_bu_b,\mathcal{E}_b)^{\mathsf{T}}$"
        "\n" "→ 边界数值通量", "#E3F2FD", BLUE)

    arrow(ax, (2.85, 4.3), (3.3, 4.3), GRAY, lw=1.8)
    arrow(ax, (2.85, 2.05), (3.3, 2.05), GRAY, lw=1.8)
    arrow(ax, (5.98, 4.3), (6.4, 4.3), GRAY, lw=1.8)
    arrow(ax, (4.65, 3.7), (4.65, 2.6), GRAY, lw=1.8)
    arrow(ax, (5.98, 2.05), (6.4, 2.05), GRAY, lw=1.8)
    arrow(ax, (8.1, 3.7), (8.1, 2.6), GRAY, lw=1.8)
    arrow(ax, (9.15, 2.05), (9.6, 2.9), GRAY, lw=1.8)

    ax.text(6.2, 0.85,
            "三条特征信息：$u_n-c$ 传入给 $p_b$；$u_n$ 与 $u_n+c$ 由内部传出\n"
            "数值示例：$u_{\\mathrm{in}}=10$ m/s， $c_{\\mathrm{in}}=300$ m/s，"
            " $c_b=299$ m/s → $J_+=1510$ m/s， $u_b=15$ m/s",
            fontsize=10.8, ha="center", color="#22303A",
            bbox=dict(boxstyle="round,pad=0.4", fc="#FAFAFA",
                      ec="#C9D2DA"))

    fig.suptitle("亚声速出口：给静压，其余由内部传出的信息补齐",
                 fontsize=15, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "fig27-subsonic-outlet-construction.png")


# ============================================ 图 28 面通量与两侧收支
def fig28():
    fig = plt.figure(figsize=(11.2, 5.0))
    ax = fig.add_axes([0.01, 0.02, 0.98, 0.845])
    ax.set_xlim(0, 11.2)
    ax.set_ylim(0, 4.35)
    ax.axis("off")

    ax.add_patch(Rectangle((1.2, 1.2), 3.6, 2.7, fc="#E3F2FD",
                           ec=BLUE, lw=1.6))
    ax.add_patch(Rectangle((6.4, 1.2), 3.6, 2.7, fc="#FFF3E0",
                           ec="#EF6C00", lw=1.6))
    ax.text(3.0, 3.0, "控制体 $L$", fontsize=14, ha="center")
    ax.text(3.0, 2.25, r"$\bar{\mathbf{U}}_L$", fontsize=14, ha="center")
    ax.text(8.2, 3.0, "控制体 $R$", fontsize=14, ha="center")
    ax.text(8.2, 2.25, r"$\bar{\mathbf{U}}_R$", fontsize=14, ha="center")

    ax.plot([5.6, 5.6], [1.05, 4.18], color="#22303A", lw=2.0)
    ax.annotate("", xy=(6.25, 4.06), xytext=(5.6, 4.06),
                arrowprops=dict(arrowstyle="-|>", color="#22303A",
                                lw=1.4, mutation_scale=11))
    ax.text(6.3, 4.16, r"$\mathbf{n}_f$（$L$ 的外法向）", fontsize=11,
            color="#22303A", va="center")

    arrow(ax, (5.1, 2.5), (6.15, 2.5), RED, lw=3.0, ms=16)
    ax.text(5.62, 2.72, r"$\widehat{\mathbf{F}}_{n,f}$", fontsize=14,
            color=RED, ha="center")
    ax.text(5.6, 0.98, r"面积 $S_f$", fontsize=11.5, ha="center",
            color=GRAY)

    ax.text(3.0, 0.55, r"$-\widehat{\mathbf{F}}_{n,f}S_f$（流出）",
            fontsize=12, ha="center", color=BLUE)
    ax.text(8.2, 0.55, r"$+\widehat{\mathbf{F}}_{n,f}S_f$（流入）",
            fontsize=12, ha="center", color="#EF6C00")

    ax.text(5.6, 0.12,
            "面法向从 $L$ 指向 $R$：对 $L$ 是流出，对 $R$ 是流入",
            fontsize=10.8, ha="center", color=GRAY)

    fig.text(0.5, 0.968, "一个面、一个通量、两个控制体",
             fontsize=15, ha="center", va="center")
    fig.text(0.5, 0.9,
             "同一份共享通量，两侧符号相反；内部面求和时相互抵消",
             fontsize=12, ha="center", va="center", color="#22303A")
    save(fig, "fig28-face-flux-two-cells.png")


# ============================================ 图 29 Rusanov 共享通量
def fig29():
    fig = plt.figure(figsize=(11.2, 5.4))
    ax = fig.add_axes([0.005, 0.01, 0.99, 0.86])
    ax.set_xlim(0, 11.2)
    ax.set_ylim(0, 5.0)
    ax.axis("off")

    def box(x, y, w, h, text, fc, ec, fs=11.5):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0.06", fc=fc,
                                    ec=ec, lw=1.5))
        ax.text(x + w / 2, y + h / 2, text, fontsize=fs,
                ha="center", va="center", color="#22303A")

    box(0.30, 3.45, 2.35, 1.15,
        "左侧状态\n" r"$\mathbf{U}_L$", "#E3F2FD", BLUE)
    box(8.55, 3.45, 2.35, 1.15,
        "右侧状态\n" r"$\mathbf{U}_R$", "#FFF3E0", "#EF6C00")
    box(3.35, 3.45, 1.95, 1.15,
        "物理通量\n" r"$\mathbf{F}(\mathbf{U}_L)$", "#E3F2FD", BLUE)
    box(5.90, 3.45, 1.95, 1.15,
        "物理通量\n" r"$\mathbf{F}(\mathbf{U}_R)$", "#FFF3E0", "#EF6C00")
    box(3.35, 2.15, 4.50, 0.90,
        r"中心平均 $\frac{1}{2}[\mathbf{F}_L+\mathbf{F}_R]$",
        "#F5F5F5", GRAY)
    box(3.35, 0.95, 4.50, 0.90,
        r"耗散修正 $-\frac{1}{2} s_{\max}(\mathbf{U}_R-\mathbf{U}_L)$",
        "#FFEBEE", RED)

    ax.text(5.60, 4.88,
            r"$s_{\max}=\max(|u_L|+c_L,\ |u_R|+c_R)$",
            fontsize=13, ha="center", va="center", color="#22303A")

    arrow(ax, (2.70, 4.02), (3.30, 4.02), GRAY, lw=1.8, ms=12)
    arrow(ax, (8.50, 4.02), (7.90, 4.02), GRAY, lw=1.8, ms=12)
    arrow(ax, (4.32, 3.40), (4.32, 3.08), GRAY, lw=1.8, ms=12)
    arrow(ax, (6.88, 3.40), (6.88, 3.08), GRAY, lw=1.8, ms=12)
    arrow(ax, (5.60, 2.10), (5.60, 1.88), GRAY, lw=1.8, ms=12)
    arrow(ax, (5.60, 0.90), (5.60, 0.44), GRAY, lw=1.8, ms=12)

    ax.text(5.60, 0.16, r"$\widehat{\mathbf{F}}=\frac{1}{2}(\mathbf{F}_L+\mathbf{F}_R)"
            r"-\frac{1}{2} s_{\max}(\mathbf{U}_R-\mathbf{U}_L)$",
            fontsize=17, ha="center", va="center", color="#22303A")
    ax.text(0.30, 0.62, "中心平均：精度", fontsize=10.5, color=GRAY)
    ax.text(0.30, 0.24, "耗散项：稳定", fontsize=10.5, color=RED)

    fig.text(0.5, 0.965,
             "Rusanov 通量：两侧信息各取一半，再加一份按波速的耗散",
             fontsize=14.8, ha="center", va="center")
    save(fig, "fig29-rusanov-flux.png")


# ============================================ 图 30 显式与隐式
def fig30():
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0))

    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.set_title("(a) 显式：用已知状态直接推进", fontsize=13)

    def bx(x, y, w, h, text, fc, ec, fs=11.5):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0.08", fc=fc,
                                    ec=ec, lw=1.5))
        ax.text(x + w / 2, y + h / 2, text, fontsize=fs,
                ha="center", va="center", color="#22303A")

    bx(0.4, 3.0, 2.5, 1.0, r"$\bar{\mathbf{U}}^{\,n}$", "#E3F2FD", BLUE)
    bx(3.7, 3.0, 2.8, 1.0, r"$L_h(\bar{\mathbf{U}}^{\,n})$", "#F5F5F5", GRAY)
    bx(7.0, 3.0, 2.5, 1.0, r"$\bar{\mathbf{U}}^{\,n+1}$", "#E8F5E9", GREEN)
    arrow(ax, (2.95, 3.5), (3.65, 3.5), GRAY, lw=1.8)
    arrow(ax, (6.55, 3.5), (6.95, 3.5), GRAY, lw=1.8)
    ax.text(5.1, 4.35, r"$\bar{\mathbf{U}}^{\,n+1}=\bar{\mathbf{U}}^{\,n}"
            r"+\Delta t\,L_h(\bar{\mathbf{U}}^{\,n})$",
            fontsize=13.5, ha="center")
    ax.text(5.0, 2.25, "每个面只算一次；\n时间步受 CFL 限制",
            fontsize=10.8, ha="center", color=GRAY)

    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.set_title("(b) 隐式：新状态也参与，需迭代", fontsize=13)

    bx(0.5, 3.0, 3.4, 1.0,
       "候选 $\\bar{\\mathbf{U}}^{\\,n+1,(k)}$", "#FFF3E0", "#EF6C00")
    bx(0.5, 1.3, 3.4, 1.0, "残差\n" r"$\mathbf{R}=\Delta\mathbf{U}-\Delta t\,L_h$",
       "#FFEBEE", RED)
    bx(5.6, 1.3, 3.6, 1.0, "线性子问题\n" r"$A\,\delta\mathbf{U}=-\mathbf{R}$",
       "#F3E5F5", "#6A1B9A")
    bx(5.6, 3.0, 3.6, 1.0, "更新\n" r"$\mathbf{U}^{(k+1)}=\mathbf{U}^{(k)}+\delta\mathbf{U}$",
       "#E8F5E9", GREEN)

    arrow(ax, (2.2, 2.95), (2.2, 2.35), GRAY, lw=1.8)
    arrow(ax, (3.95, 1.8), (5.55, 1.8), GRAY, lw=1.8)
    arrow(ax, (7.4, 2.35), (7.4, 2.95), GRAY, lw=1.8)
    arrow(ax, (5.55, 3.5), (3.95, 3.5), GRAY, lw=1.8)
    ax.text(4.75, 3.65, "迭代", fontsize=10.5, ha="center", color=GRAY)
    ax.text(5.0, 0.55, r"收敛后保存为 $\bar{\mathbf{U}}^{\,n+1}$；"
            "同一物理步内迭代", fontsize=10.8, ha="center", color=GRAY)

    fig.suptitle("同一个空间算子 $L_h$，两种时间组织方式",
                 fontsize=15, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "fig30-explicit-vs-implicit.png")


# ============================================ 图 31 低马赫数时间尺度
def fig31():
    labels = [r"显式步长 $\Delta t$（$h=1$ cm）",
              r"声学穿越 $L/c$（$L=1$ km）",
              r"流动穿越 $L/u$（$L=1$ km）"]
    vals = [8.6e-6, 2.0, 67.0]
    cols = [RED, "#EF6C00", GREEN]

    fig, ax = plt.subplots(figsize=(11.0, 4.6))
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, color=cols, height=0.5, alpha=0.9)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=11.5)
    ax.set_xlabel("时间尺度 [s]（对数轴）", fontsize=12)
    for yi, v in zip(y, vals):
        if v < 1e-3:
            expo = int(np.floor(np.log10(v)))
            mant = v / 10.0 ** expo
            txt = r"$%.1f\times10^{%d}$ s" % (mant, expo)
        else:
            txt = "%.0f s" % v
        ax.text(v * 1.25, yi, txt, va="center", fontsize=11)
    ax.set_xlim(1e-6, 300)
    ax.grid(True, axis="x", which="both", color="#EEF2F5", lw=0.8)
    ax.set_axisbelow(True)
    ax.text(0.97, 0.96,
            "甲烷：30 MPa、380 K， $c\\approx508$ m/s； $u=15$ m/s，"
            " $\\mathrm{Ma}\\approx0.03$\n"
            "声学步长比流动时间小 4 个量级 → 显式可压缩在低马赫数下很慢",
            transform=ax.transAxes, fontsize=10.8, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec="#C9D2DA"))
    fig.suptitle("低马赫数：慢流动被快声波拖住", fontsize=15, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "fig31-low-mach-timescales.png")


# ============================================ 图 32 Sod 数值 vs 精确解
def fig32():
    import sys
    sys.path.insert(0, HERE)
    import fluid_lib as fl

    t = 0.2
    runs = [(200, "#1E88E5", 1.4), (800, "#C62828", 1.4)]
    data = []
    for n, col, lw in runs:
        r = fl.fv_sod(n, t, cfl=0.45)
        U = r["U"]
        prim = np.array([fl.primitive(u) for u in U])
        data.append((r["x"], prim, col, n))

    xf = np.linspace(0, 1, 1200)
    rex, uex, pex = fl.sod_exact(xf, t)

    fig, axes = plt.subplots(3, 1, figsize=(10.6, 8.2), sharex=True)
    fields = [(0, r"密度 $\rho$"), (1, r"速度 $u$"), (2, r"压力 $p$")]
    exact = [rex, uex, pex]
    for ax, (k, lab), ex in zip(axes, fields, exact):
        ax.plot(xf, ex, "-", color="#22303A", lw=2.0,
                label="精确解")
        for x, prim, col, n in data:
            ax.plot(x, prim[:, k], "-", color=col, lw=1.3,
                    label="%d 个单元" % n)
        ax.set_ylabel(lab, fontsize=12)
        ax.grid(True, color="#EEF2F5", lw=0.8)
        ax.set_axisbelow(True)
        ax.legend(loc="best", fontsize=9.5)

    ps, us, cL, cR = fl.sod_star_state()
    cstar = cL * (ps / 1.0) ** (0.4 / 2.8)
    S = cR * np.sqrt(2.4 / 2.8 * ps / 0.1 + 0.4 / 2.8)
    axes[0].axvline(0.5 + (us - cstar) * t, color=GRAY, ls=":", lw=1.0)
    axes[0].axvline(0.5 + us * t, color=GRAY, ls=":", lw=1.0)
    axes[0].axvline(0.5 + S * t, color=GRAY, ls=":", lw=1.0)
    axes[0].text(0.42, 0.88, "稀疏波", fontsize=9.5, color=GRAY)
    axes[0].text(0.64, 0.88, "接触", fontsize=9.5, color=GRAY)
    axes[0].text(0.79, 0.88, "激波", fontsize=9.5, color=GRAY)
    axes[0].set_ylim(0, 1.05)
    axes[2].set_xlabel("位置 $x$", fontsize=12)
    fig.suptitle("Sod 激波管：一阶 Rusanov 也能抓住三道结构，代价是抹宽",
                 fontsize=14.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, "fig32-sod-shocktube.png")


# ============================================ 图 33 激波跃变控制体
def fig33():
    fig = plt.figure(figsize=(11.2, 4.9))
    ax = fig.add_axes([0.005, 0.01, 0.99, 0.855])
    ax.set_xlim(0, 11.2)
    ax.set_ylim(0, 4.7)
    ax.axis("off")

    ax.add_patch(Rectangle((0.4, 1.5), 5.1, 2.6, fc="#E3F2FD",
                           ec=BLUE, lw=1.4))
    ax.add_patch(Rectangle((5.5, 1.5), 5.3, 2.6, fc="#FFF3E0",
                           ec="#EF6C00", lw=1.4))
    ax.text(2.6, 3.35, "波前 1\n" r"$\rho_1,\ u_1,\ p_1$", fontsize=13,
            ha="center", va="center")
    ax.text(8.4, 3.35, "波后 2\n" r"$\rho_2,\ u_2,\ p_2$", fontsize=13,
            ha="center", va="center")

    ax.add_patch(Rectangle((5.15, 1.15), 0.5, 3.1, fc="#FFEBEE",
                           ec=RED, lw=2.0))
    ax.annotate("", xy=(7.2, 4.5), xytext=(5.4, 4.5),
                arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.2,
                                mutation_scale=14))
    ax.text(7.35, 4.5, r"激波速度 $S$（固定坐标）", fontsize=12,
            color=RED, va="center")

    arrow(ax, (3.9, 1.95), (5.1, 1.95), GREEN, lw=2.2, ms=12)
    ax.text(4.5, 2.14, r"相对流入 $S-u_1$", fontsize=10.8,
            color=GREEN, ha="center")
    arrow(ax, (5.75, 1.95), (6.95, 1.95), "#6A1B9A", lw=2.2, ms=12)
    ax.text(6.35, 2.14, r"相对流出 $S-u_2$", fontsize=10.8,
            color="#6A1B9A", ha="center")

    ax.text(5.6, 0.62,
            r"$\rho_1(S-u_1)=\rho_2(S-u_2)$"
            r"$\qquad$ 每单位面积的质量通过率相同",
            fontsize=14.5, ha="center", color="#22303A")
    ax.text(5.6, 0.18,
            "动量与能量同理：把压力、内能与动能一并计入，得到 Rankine–Hugoniot 条件",
            fontsize=10.5, ha="center", color=GRAY)

    fig.text(0.5, 0.965,
             "随激波移动的控制体：状态可以跳，收支不能断",
             fontsize=15, ha="center", va="center")
    save(fig, "fig33-shock-jump-control-volume.png")


# ============================================================ 图 34 单元家族
def fig34():
    """三类单元与 Abaqus 命名的对照。"""
    fig, ax = plt.subplots(figsize=(11.4, 7.4))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8.5)
    ax.axis("off")

    def box(x, y, w, h, title, lines, color, fc):
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.06,rounding_size=0.12",
            fc=fc, ec=color, lw=1.5))
        ax.text(x + w / 2, y + h - 0.24, title, ha="center", va="top",
                fontsize=10.6, color=color, weight="bold")
        for k, s in enumerate(lines):
            ax.text(x + 0.18, y + h - 0.60 - 0.285 * k, s,
                    ha="left", va="top", fontsize=9.0, color="#22303A")

    ax.add_patch(FancyBboxPatch((2.4, 7.35), 7.2, 0.88,
                                boxstyle="round,pad=0.08,rounding_size=0.16",
                                fc="#E8F1FB", ec=BLUE, lw=1.8))
    ax.text(6.0, 7.98, "单元类型 = 未知量的选择", ha="center", va="center",
            fontsize=14.5, color=BLUE, weight="bold")
    ax.text(6.0, 7.55, "名字是编码：族 + 维数 + 节点数 + 积分 / 扩展",
            ha="center", va="center", fontsize=10.2, color="#22303A")

    heads = [
        (1.95, "三维实体", "3 个平动 / 节点", GREEN),
        (6.0, "结构单元", "平动 + 转动", BLUE),
        (10.05, "二维降维", "3 个平动 / 节点", RED),
    ]
    for x, t, s, c in heads:
        ax.text(x, 6.78, t, ha="center", va="center", fontsize=12.2,
                color=c, weight="bold")
        ax.text(x, 6.48, s, ha="center", va="center", fontsize=9.6, color=GRAY)
        arrow(ax, (x, 6.34), (x, 6.18), c, lw=1.4, ms=10)

    box(0.3, 4.55, 3.3, 1.72, "C3D8 / C3D8R / C3D8I",
        ["8 节点六面体，3 平动/节点",
         "全积分 / 减缩 / 非协调模式",
         "弯曲优先 C3D8I 或 C3D10"], GREEN, "#E8F5E9")
    box(0.3, 2.62, 3.3, 1.72, "C3D10 / C3D4",
        ["10 节点二次四面体（边中点）",
         "4 节点线性四面体：常应变",
         "复杂几何常用 C3D10"], GREEN, "#E8F5E9")
    box(0.3, 0.69, 3.3, 1.72, "C3D8P / C3D8T",
        ["附加孔压 / 温度自由度",
         "每节点 4 个未知量",
         "Soils 步用非对称求解"], "#6A1B9A", "#F3E5F5")

    box(4.25, 4.55, 3.5, 1.72, "S4 / S4R / S8R",
        ["4 或 8 节点壳，6 自由度",
         "厚度 t 是截面输入，不划网格",
         "薄：t 小于跨度约 1/10"], BLUE, "#E3F2FD")
    box(4.25, 2.62, 3.5, 1.72, "B31 / B32",
        ["2 / 3 节点梁，6 自由度",
         "输入 A、I、J，输出截面内力",
         "细长：截面小于长度约 1/10"], BLUE, "#E3F2FD")
    box(4.25, 0.69, 3.5, 1.72, "T3D2 / T3D3",
        ["2 / 3 节点杆，只有 3 平动",
         "只受轴力，不能受弯",
         "辐条、拉索、桁架腹杆"], "#00838F", "#E0F7FA")

    box(8.2, 4.55, 3.5, 1.72, "CPS4 / CPS8",
        [r"平面应力：$\sigma_z=0$",
         "薄片在自身平面内受力",
         "链片、膜、薄板面内"], RED, "#FFEBEE")
    box(8.2, 2.62, 3.5, 1.72, "CPE4 / CPE8",
        [r"平面应变：$\varepsilon_z=0$",
         "长隧道、长坝、长边坡",
         r"$\sigma_z=\nu(\sigma_x+\sigma_y)$"], RED, "#FFEBEE")
    box(8.2, 0.69, 3.5, 1.72, "CAX4 / CAX8",
        [r"轴对称：$\varepsilon_\theta=u_r/r$",
         "只画半个剖面",
         "套管、储罐、井筒"], RED, "#FFEBEE")

    ax.add_patch(FancyBboxPatch((0.3, 0.06), 11.4, 0.46,
                                boxstyle="round,pad=0.04,rounding_size=0.08",
                                fc="#F5F5F5", ec=GRAY, lw=1.0))
    ax.text(6.0, 0.29,
            "C3D8R = Continuum + 3D + 8 nodes + Reduced    |    "
            "S4R = Shell + 4 nodes + Reduced    |    "
            "C3D10 = Continuum + 3D + 10 nodes",
            ha="center", va="center", fontsize=9.4, color="#22303A")

    fig.text(0.5, 0.975, "单元家族：名字里写着自由度与适用比例",
             fontsize=15.5, ha="center", va="center")
    save(fig, "fig34-element-family.png")


# ============================================================ 图 35 四类单元
def fig35():
    """实体、壳、梁、杆的几何与自由度。"""
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.8))
    for ax in axes.ravel():
        ax.set_aspect("equal")
        ax.axis("off")

    # (a) 实体
    ax = axes[0, 0]
    ax.set_xlim(-0.1, 2.3)
    ax.set_ylim(-0.62, 2.2)
    d = np.array([0.42, 0.42])
    front = np.array([[0.2, 0.2], [1.45, 0.2], [1.45, 1.45], [0.2, 1.45]])
    back = front + d
    ax.add_patch(Polygon(back, closed=True, fc="#E8F5E9", ec=GREEN, lw=1.6))
    for p, q in zip(front, back):
        ax.plot([p[0], q[0]], [p[1], q[1]], color=GREEN, lw=1.0)
    ax.add_patch(Polygon(front, closed=True, fc="#C8E6C9", ec=GREEN, lw=2.0))
    ax.scatter(np.r_[front[:, 0], back[:, 0]], np.r_[front[:, 1], back[:, 1]],
               s=34, color=GREEN, zorder=5)
    ax.text(0.80, 0.62, "8 个节点", ha="center", fontsize=9.5, color=GREEN)
    arrow(ax, (1.45, 1.45), (2.0, 1.45), BLUE, lw=1.8, ms=11)
    arrow(ax, (1.45, 1.45), (1.45, 2.0), GREEN, lw=1.8, ms=11)
    arrow(ax, (1.45, 1.45), (1.82, 1.82), RED, lw=1.8, ms=11)
    ax.text(2.02, 1.40, r"$u_x$", fontsize=10, color=BLUE)
    ax.text(1.43, 2.02, r"$u_y$", fontsize=10, color=GREEN)
    ax.text(1.78, 1.86, r"$u_z$", fontsize=10, color=RED)
    ax.set_title("(a) 实体：C3D8 / C3D10", fontsize=12, color=GREEN)
    ax.text(-0.05, -0.58, "3 个平动 / 节点；应力在单元内积分",
            fontsize=9.4, color="#22303A", va="bottom")

    # (b) 壳
    ax = axes[0, 1]
    ax.set_xlim(-0.12, 1.75)
    ax.set_ylim(-0.2, 1.25)
    x = np.linspace(0.08, 1.55, 120)
    y = 0.22 * np.sin(np.pi * (x - 0.08) / 1.47)
    y0 = 0.62
    top = y0 + y + 0.07
    bot = y0 + y - 0.07
    ax.fill_between(x, bot, top, color="#BBDEFB", ec=BLUE, lw=1.7)
    ax.plot(x, y0 + y, ls="--", color=BLUE, lw=1.2)
    for xx in np.linspace(0.18, 1.45, 5):
        yy = y0 + 0.22 * np.sin(np.pi * (xx - 0.08) / 1.47)
        ax.scatter([xx], [yy], s=30, color=BLUE, zorder=5)
        arrow(ax, (xx, yy), (xx, yy + 0.14), RED, lw=1.1, ms=8)
    ax.annotate("", xy=(1.22, y0 + 0.30), xytext=(1.22, y0 + 0.16),
                arrowprops=dict(arrowstyle="<->", color=GREEN, lw=1.2))
    ax.text(1.26, y0 + 0.20, r"$t$", fontsize=10, color=GREEN)
    ax.text(0.78, 0.30, "中面铺网格，厚度是输入", ha="center",
            fontsize=9.6, color="#22303A")
    ax.set_title("(b) 壳：S4 / S4R / S8R", fontsize=12, color=BLUE)

    # (c) 梁
    ax = axes[1, 0]
    ax.set_xlim(-0.1, 2.1)
    ax.set_ylim(-0.1, 1.45)
    ax.plot([0.15, 1.55], [0.62, 0.62], color="#1565C0", lw=3.2)
    for xx in np.linspace(0.15, 1.55, 5):
        ax.scatter([xx], [0.62], s=34, color=BLUE, zorder=5)
    ax.add_patch(Ellipse((1.55, 0.62), 0.18, 0.62, fc="#BBDEFB",
                         ec=BLUE, lw=1.5))
    ax.annotate("", xy=(1.55, 0.93), xytext=(1.55, 0.31),
                arrowprops=dict(arrowstyle="<->", color=GREEN, lw=1.1))
    ax.text(1.72, 0.58, "D", fontsize=10, color=GREEN)
    ax.add_patch(Arc((1.2, 0.62), 0.32, 0.32, theta1=200, theta2=340,
                     color=RED, lw=1.8))
    arrow(ax, (1.35, 0.53), (1.37, 0.47), RED, lw=1.4, ms=9)
    ax.text(1.15, 0.34, r"$\theta$", fontsize=10, color=RED)
    ax.text(0.82, 0.08, "一根线 + 截面性质 A、I、J", ha="center",
            fontsize=9.6, color="#22303A")
    ax.set_title("(c) 梁：B31 / B32", fontsize=12, color=BLUE)

    # (d) 杆
    ax = axes[1, 1]
    ax.set_xlim(-0.1, 2.1)
    ax.set_ylim(-0.1, 1.45)
    ax.plot([0.25, 1.55], [0.7, 0.7], color="#00838F", lw=3.2)
    ax.scatter([0.25, 1.55], [0.7, 0.7], s=40, color="#00838F", zorder=5)
    arrow(ax, (0.25, 0.7), (-0.02, 0.7), GREEN, lw=1.8, ms=11)
    arrow(ax, (1.55, 0.7), (1.85, 0.7), GREEN, lw=1.8, ms=11)
    ax.text(0.55, 0.80, r"$N$", fontsize=12, color=GREEN)
    ax.plot([1.18, 1.33], [0.94, 1.09], color=RED, lw=2.2)
    ax.plot([1.18, 1.33], [1.09, 0.94], color=RED, lw=2.2)
    ax.text(1.37, 1.01, "无弯曲", fontsize=9.5, color=RED, va="center")
    ax.text(0.9, 0.12, "只有轴力，节点不转动", ha="center",
            fontsize=9.6, color="#22303A")
    ax.set_title("(d) 杆：T3D2 / T3D3", fontsize=12, color="#00838F")

    fig.suptitle("四类单元：几何决定未知量", fontsize=15.5, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, "fig35-four-element-types.png")


# ============================================================ 图 36 悬臂梁
def fig36():
    """悬臂梁：Q4 / Q8 / T3 的挠曲线与端部误差。"""
    exact = el.beam_tip_exact()
    cases = [
        ("Q4 10×1", el.mesh_q4(10, 1),
         lambda c: el.isoparametric_stiffness(c, el.q4_shape, el.gauss_2d(2)),
         RED, "-"),
        ("Q8 10×1", el.mesh_q8(10, 1),
         lambda c: el.isoparametric_stiffness(c, el.q8_shape, el.gauss_2d(3)),
         BLUE, "-"),
        ("T3 10×1", el.mesh_t3(10, 1), el.cst_stiffness, GREEN, "-"),
        ("Q4 80×8（640 块）", el.mesh_q4(80, 8),
         lambda c: el.isoparametric_stiffness(c, el.q4_shape, el.gauss_2d(2)),
         "#EF6C00", "--"),
    ]

    fig = plt.figure(figsize=(8.9, 8.7))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.45, 1.0], hspace=0.34)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    xe, de = el.beam_centerline_exact()
    ax1.plot(xe / el.BEAM_L, de / exact, color="black", lw=2.8,
             label="Euler–Bernoulli 真值", zorder=1)

    ratios = []
    for label, (nodes, elems), fn, color, ls in cases:
        U = el.assemble_and_solve(nodes, elems, fn)
        x, u = el.centerline(nodes, U)
        ratio = u[-1] / exact
        ratios.append((label, ratio))
        ax1.plot(x / el.BEAM_L, u / exact, color=color, ls=ls, lw=2.2,
                 label="%s  %.3f" % (label, ratio), zorder=3)

    ax1.axhline(1.0, color=GRAY, ls=":", lw=1.0)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1.06)
    ax1.set_xlabel(r"$x/L$", fontsize=11)
    ax1.set_ylabel("挠度 / 梁理论真值", fontsize=11)
    ax1.set_title("同一根梁、同样 10 块：单元决定误差", fontsize=13.5)
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=9.2, loc="upper left", framealpha=0.95)
    ax1.text(0.02, 0.03,
             "P = 1000 N，L = 100 mm，H = 10 mm，t = 10 mm；真值 1.9048 mm",
             fontsize=8.8, color=GRAY, transform=ax1.transAxes)

    labels = ["Q4\n10×1", "Q4\n10×2", "Q4\n10×4", "Q4\n10×8",
              "Q4\n80×8", "Q8\n10×1", "T3\n10×1"]
    vals = [ratios[0][1], 0.71109, 0.72015, 0.72247,
            ratios[3][1], ratios[1][1], ratios[2][1]]
    colors = [RED, RED, RED, RED, "#EF6C00", BLUE, GREEN]
    bars = ax2.bar(range(len(vals)), vals, color=colors, width=0.66)
    ax2.axhline(1.0, color="black", ls="--", lw=1.4)
    ax2.set_xticks(range(len(vals)))
    ax2.set_xticklabels(labels, fontsize=9.2)
    ax2.set_ylim(0, 1.12)
    ax2.set_ylabel("端部挠度 / 真值", fontsize=11)
    ax2.set_title("端部精度：全积分、二次与线性三角形", fontsize=12.5)
    ax2.grid(axis="y", alpha=0.3)
    for b, v in zip(bars, vals):
        ax2.text(b.get_x() + b.get_width() / 2, v + 0.025, "%.1f%%" % (100 * v),
                 ha="center", fontsize=9.0, color="#22303A")
    ax2.text(0.01, 0.92, "虚线 = 100%", fontsize=9.0, color=GRAY,
             transform=ax2.transAxes)

    fig.text(0.5, 0.985, "悬臂梁基准：为什么厚度方向的层数不是随便定的",
             fontsize=15.0, ha="center", va="center")
    save(fig, "fig36-cantilever-accuracy.png")


# ============================================================ 图 37 剪切锁死
def fig37():
    """纯弯曲下 Q4 的寄生剪切。"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.8),
                                   gridspec_kw={"width_ratios": [1.05, 1.0]})

    # 左：单元变形
    a, h, kap = 1.25, 0.30, 1.4
    ax1.set_aspect("equal")
    ax1.set_xlim(-0.18, 1.55)
    ax1.set_ylim(-0.14, 0.72)
    ax1.axis("off")
    ax1.add_patch(Rectangle((0, 0), a, h, fill=False, ls="--",
                            ec=GRAY, lw=1.3))
    # 精确弯曲：上下边是弧（中性轴在中面 y=h/2）
    xs = np.linspace(0, a, 80)
    for yy in [0.0, h]:
        ax1.plot(xs + 0.10 * (-kap * xs * (yy - h / 2)),
                 yy + 0.10 * (0.5 * kap * xs ** 2), color=GREEN,
                 lw=1.4, ls=":", zorder=2)
    # 双线性变形：顶点用同样的节点位移
    corners = np.array([[0, 0], [a, 0], [a, h], [0, h]], float)
    ux = -kap * corners[:, 0] * (corners[:, 1] - h / 2)
    uy = 0.5 * kap * corners[:, 0] ** 2
    defo = corners + 0.10 * np.c_[ux, uy]
    ax1.add_patch(Polygon(defo, closed=True, fc="#BBDEFB", ec=BLUE,
                          lw=2.0, alpha=0.75, zorder=3))
    for p, q in zip(corners, defo):
        ax1.plot([p[0], q[0]], [p[1], q[1]], color=BLUE, lw=0.8, ls="-.")
    ax1.scatter(corners[:, 0], corners[:, 1], s=26, color=BLUE, zorder=6)
    g = 1 / np.sqrt(3)
    for xi in [-g, g]:
        for eta in [-g, g]:
            xg = a / 2 * (1 + xi)
            yg = h / 2 * (1 + eta)
            ax1.scatter([xg], [yg], s=42, marker="s", color=RED, zorder=7)
            dy = 0.055 if xi * eta > 0 else -0.055
            arrow(ax1, (xg, yg), (xg, yg + dy), RED, lw=1.2, ms=7)
    ax1.text(a / 2, h + 0.10, "直边只能变梯形 → 寄生的 $\\gamma_{xy}$",
             ha="center", fontsize=10.2, color=BLUE)
    ax1.text(-0.14, 0.05, "虚线：精确弯曲", fontsize=8.8, color=GREEN,
             rotation=90, va="bottom")
    ax1.text(a + 0.05, -0.10, "红方块：2×2 Gauss 点", fontsize=8.8,
             color=RED)
    ax1.set_title("(a) 一个 Q4 单元的纯弯曲", fontsize=12.5)

    # 右：寄生剪切随长厚比
    ar = np.logspace(0, 2, 200)
    ax2.loglog(ar, ar / (2 * np.sqrt(3)), color=BLUE, lw=2.4,
               label=r"解析：$\gamma=\kappa a/(2\sqrt{3})$")
    for aa in [1, 10, 100]:
        ax2.loglog([aa], [el.parasitic_shear(1.0, aa)], "o", color=RED,
                   ms=7, zorder=5)
    ax2.axhline(1.0, color=GRAY, ls=":", lw=1.0)
    ax2.text(1.15, 1.12, r"参考线 $\gamma=\kappa h$", fontsize=9.0, color=GRAY)
    ax2.annotate(r"$a/h=10$：$\gamma=2.89\,\kappa h$",
                 xy=(10, 2.89), xytext=(1.8, 9.0), fontsize=9.6,
                 color=RED,
                 arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.2))
    ax2.set_xlabel("长厚比 $a/h$（单元长度 / 厚度）", fontsize=10.5)
    ax2.set_ylabel(r"寄生剪切 $\gamma/(\kappa h)$", fontsize=10.5)
    ax2.set_title("(b) 长厚比越大，锁死越重", fontsize=12.5)
    ax2.grid(which="both", alpha=0.3)
    ax2.legend(fontsize=9.2, loc="upper left")

    fig.suptitle("剪切锁死：直边弯不了，只能靠假剪切消化弯曲",
                 fontsize=14.5, y=1.0)
    fig.tight_layout()
    save(fig, "fig37-shear-locking.png")


# ============================================================ 图 38 沙漏
def fig38():
    """减缩积分的零能模式与特征值谱。"""
    efull = el.element_eigenvalues(el.q4_shape, el.gauss_2d(2))
    ered = el.element_eigenvalues(el.q4_shape, el.gauss_2d(1))
    modes = el.hourglass_modes()

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.4),
                             gridspec_kw={"width_ratios": [1.2, 1, 1]})
    ax = axes[0]
    floor = 1e-14
    idx = np.arange(1, 9)
    ax.plot(idx - 0.08, np.log10(np.maximum(efull / efull.max(), floor)),
            "o", color=BLUE, ms=8, label="Q4 2×2（全积分）")
    ax.plot(idx + 0.08, np.log10(np.maximum(ered / ered.max(), floor)),
            "x", color=RED, ms=9, mew=2.2, label="Q4 1×1（减缩）")
    ax.axhline(-10, color=GRAY, ls=":", lw=1.0)
    ax.text(0.55, -9.6, "数值零", fontsize=8.8, color=GRAY)
    ax.annotate("3 个刚体模态", xy=(2, -12.5), xytext=(2.7, -6.0),
                fontsize=9.2, color=BLUE,
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.1))
    ax.annotate("2 个沙漏模态", xy=(4.6, -12.5), xytext=(5.0, -5.2),
                fontsize=9.2, color=RED,
                arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.1))
    ax.set_xlim(0.4, 8.6)
    ax.set_ylim(-14.5, 0.6)
    ax.set_xlabel("模态序号", fontsize=10.5)
    ax.set_ylabel(r"$\log_{10}(\lambda/\lambda_{\max})$", fontsize=10.5)
    ax.set_title("(a) 单刚特征值谱", fontsize=12)
    ax.legend(fontsize=8.8, loc="lower left")
    ax.grid(alpha=0.25)

    xy = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    for k, (ax, mode, color, name) in enumerate([
            (axes[1], modes[0], BLUE, "沙漏模态 1：x 方向交替"),
            (axes[2], modes[1], RED, "沙漏模态 2：y 方向交替")]):
        ax.set_aspect("equal")
        ax.set_xlim(-0.35, 1.35)
        ax.set_ylim(-0.35, 1.35)
        ax.axis("off")
        ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, ls="--",
                               ec=GRAY, lw=1.2))
        ux = mode[0::2]
        uy = mode[1::2]
        scale = 0.30
        ax.quiver(xy[:, 0], xy[:, 1], ux, uy, angles="xy",
                  scale_units="xy", scale=1 / scale, color=color,
                  width=0.012, zorder=5)
        # 双线性变形后的边界
        ts = np.linspace(0, 1, 40)
        edges = [np.c_[ts, np.zeros_like(ts)],
                 np.c_[np.ones_like(ts), ts],
                 np.c_[1 - ts, np.ones_like(ts)],
                 np.c_[np.zeros_like(ts), 1 - ts]]
        for edge in edges:
            xi = 2 * edge[:, 0] - 1
            eta = 2 * edge[:, 1] - 1
            N, _ = el.q4_shape(xi, eta)
            dx = (scale * ux) @ N
            dy = (scale * uy) @ N
            ax.plot(edge[:, 0] + dx, edge[:, 1] + dy, color=color, lw=2.2)
        ax.scatter(xy[:, 0], xy[:, 1], s=30, color=color, zorder=6)
        ax.set_title("(%s) %s" % ("bc"[k], name), fontsize=11.2, color=color)

    fig.suptitle("减缩积分：消掉了寄生剪切，也放进了零能模式",
                 fontsize=14.5, y=1.02)
    fig.tight_layout()
    save(fig, "fig38-hourglass-modes.png")


# ============================================================ 图 39 二维降维
def fig39():
    """平面应力、平面应变与轴对称的运动学约束。"""
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.5))
    for ax in axes:
        ax.set_aspect("equal")
        ax.axis("off")

    # (a) 平面应力
    ax = axes[0]
    ax.set_xlim(-0.15, 1.65)
    ax.set_ylim(-0.15, 1.52)
    d = np.array([0.18, 0.18])
    front = np.array([[0.1, 0.1], [1.2, 0.1], [1.2, 0.95], [0.1, 0.95]])
    back = front + d
    ax.add_patch(Polygon(back, closed=True, fc="#FFCDD2", ec=RED, lw=1.3))
    for p, q in zip(front, back):
        ax.plot([p[0], q[0]], [p[1], q[1]], color=RED, lw=0.9)
    ax.add_patch(Polygon(front, closed=True, fc="#FFEBEE", ec=RED, lw=1.8))
    arrow(ax, (0.55, 0.52), (1.0, 0.52), GREEN, lw=1.6, ms=10)
    arrow(ax, (1.0, 0.52), (0.55, 0.52), GREEN, lw=1.6, ms=10)
    arrow(ax, (0.62, 0.30), (0.62, 0.78), BLUE, lw=1.6, ms=10)
    arrow(ax, (0.62, 0.78), (0.62, 0.30), BLUE, lw=1.6, ms=10)
    ax.annotate("", xy=(1.38, 0.28), xytext=(1.38, 0.10),
                arrowprops=dict(arrowstyle="<->", color=GRAY, lw=1.0))
    ax.text(1.42, 0.17, r"$t$", fontsize=10, color=GRAY)
    ax.text(0.72, -0.05, r"$\sigma_z=0$", fontsize=11, color=RED,
            ha="center")
    ax.text(0.72, 1.26, r"$\varepsilon_z=-\nu(\sigma_x+\sigma_y)/E$",
            fontsize=9.8, ha="center", color="#22303A")
    ax.set_title("(a) 平面应力：薄片", fontsize=12.5, color=RED)

    # (b) 平面应变
    ax = axes[1]
    ax.set_xlim(-0.15, 1.75)
    ax.set_ylim(-0.15, 1.62)
    d = np.array([0.42, 0.22])
    front = np.array([[0.1, 0.18], [1.0, 0.18], [1.0, 1.0], [0.1, 1.0]])
    back = front + d
    ax.add_patch(Polygon(back, closed=True, fc="#FFCDD2", ec=RED, lw=1.3))
    for p, q in zip(front, back):
        ax.plot([p[0], q[0]], [p[1], q[1]], color=RED, lw=0.9)
    ax.add_patch(Polygon(front, closed=True, fc="#FFEBEE", ec=RED, lw=1.8))
    for x0 in [0.25, 0.55, 0.85]:
        ax.plot([x0, x0 + d[0]], [1.0, 1.22], color=RED, lw=0.6, ls=":")
    ax.text(0.55, 0.50, r"$\varepsilon_z=0$", fontsize=10.5, color=RED,
            ha="center")
    ax.annotate("", xy=(1.52, 0.20), xytext=(1.52, 0.70),
                arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.2))
    ax.text(1.56, 0.42, r"$z$", fontsize=10, color=GRAY)
    ax.text(0.62, 1.40, r"$\sigma_z=\nu(\sigma_x+\sigma_y)$",
            fontsize=9.8, ha="center", color="#22303A")
    ax.set_title("(b) 平面应变：长柱体", fontsize=12.5, color=RED)

    # (c) 轴对称
    ax = axes[2]
    ax.set_xlim(-0.15, 1.5)
    ax.set_ylim(-0.1, 1.3)
    ax.plot([0, 0], [-0.05, 1.25], ls="-.", color=GRAY, lw=1.2)
    ax.text(-0.06, 1.20, "轴", fontsize=9.5, color=GRAY, ha="right")
    ax.add_patch(Rectangle((0.18, 0.12), 0.75, 0.95, fc="#FFEBEE",
                           ec=RED, lw=1.8))
    for z in [0.35, 0.65, 0.95]:
        ax.annotate("", xy=(0.18, z), xytext=(0.02, z),
                    arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.2))
    ax.text(0.08, 1.02, r"$u_r$", fontsize=10, color=GREEN)
    ax.text(0.55, 0.58, "半个剖面\n绕轴转一圈", ha="center", fontsize=9.2,
            color=RED)
    ax.text(0.55, 0.00, r"$\varepsilon_\theta=u_r/r$", fontsize=10.5,
            ha="center", color="#22303A")
    ax.set_title("(c) 轴对称：回转体", fontsize=12.5, color=RED)

    fig.suptitle("降维不是少画一个方向，而是给那个方向加约束",
                 fontsize=14.5, y=1.02)
    fig.tight_layout()
    save(fig, "fig39-2d-reductions.png")


# ============================================================ 图 40 选型流程
def fig40():
    """按几何比例与物理扩展选单元。"""
    fig, ax = plt.subplots(figsize=(10.8, 8.4))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis("off")

    def diamond(cx, cy, w, h, text):
        pts = [(cx, cy + h / 2), (cx + w / 2, cy),
               (cx, cy - h / 2), (cx - w / 2, cy)]
        ax.add_patch(Polygon(pts, closed=True, fc="#FFF8E1",
                             ec="#F9A825", lw=1.6))
        ax.text(cx, cy, text, ha="center", va="center", fontsize=9.8)

    def result(x, y, w, h, title, lines, color, fc):
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.10",
            fc=fc, ec=color, lw=1.5))
        ax.text(x + w / 2, y + h - 0.22, title, ha="center", va="top",
                fontsize=11.0, color=color, weight="bold")
        for k, s in enumerate(lines):
            ax.text(x + w / 2, y + h - 0.55 - 0.26 * k, s,
                    ha="center", va="top", fontsize=9.0, color="#22303A")

    ax.add_patch(FancyBboxPatch((3.0, 9.05), 6.0, 0.72,
                                boxstyle="round,pad=0.06,rounding_size=0.12",
                                fc="#E8F1FB", ec=BLUE, lw=1.6))
    ax.text(6.0, 9.41, "选单元：先看几何比例，再看物理扩展",
            ha="center", va="center", fontsize=12.5, color=BLUE,
            weight="bold")

    diamond(4.4, 7.85, 3.2, 1.05,
            "一个方向很薄？\n$t<L/10$")
    diamond(4.4, 5.95, 3.2, 1.05,
            "截面小于长度 1/10？")
    diamond(4.4, 4.05, 3.2, 1.05, "只传递轴力？")

    result(7.4, 7.35, 4.3, 1.00, "壳：S4 / S4R / S8R",
           ["6 自由度；厚度是截面输入", "中面网格，不必切厚度层"],
           BLUE, "#E3F2FD")
    result(7.4, 5.45, 4.3, 1.00, "梁：B31 / B32",
           ["6 自由度；输入 A、I、J", "输出截面轴力、剪力、弯矩"],
           BLUE, "#E3F2FD")
    result(7.4, 3.55, 4.3, 1.00, "杆：T3D2 / T3D3",
           ["3 平动；只受轴力", "不能受弯，不能受压失稳"],
           "#00838F", "#E0F7FA")
    result(2.4, 1.95, 7.2, 1.10, "三维实体：C3D8 / C3D8R / C3D10",
           ["弯曲优先 C3D8I、C3D10；C3D8 全积分偏硬",
            "复杂几何用 C3D10，规则六面体用 C3D8R + 足够层数"],
           GREEN, "#E8F5E9")

    arrow(ax, (6.0, 9.02), (4.4, 8.42), GRAY, lw=1.5, ms=11)
    arrow(ax, (6.0, 7.32), (7.4, 7.85), GREEN, lw=1.6, ms=11)
    ax.text(6.55, 7.76, "是", fontsize=10, color=GREEN)
    arrow(ax, (4.4, 7.32), (4.4, 6.50), GRAY, lw=1.5, ms=11)
    ax.text(4.52, 6.90, "否", fontsize=10, color=GRAY)
    arrow(ax, (6.0, 5.42), (7.4, 5.95), GREEN, lw=1.6, ms=11)
    ax.text(6.55, 5.86, "是", fontsize=10, color=GREEN)
    arrow(ax, (4.4, 5.42), (4.4, 4.60), GRAY, lw=1.5, ms=11)
    ax.text(4.52, 5.00, "否", fontsize=10, color=GRAY)
    arrow(ax, (6.0, 3.52), (7.4, 4.05), GREEN, lw=1.6, ms=11)
    ax.text(6.55, 3.96, "是", fontsize=10, color=GREEN)
    arrow(ax, (4.4, 3.52), (4.4, 3.08), GRAY, lw=1.5, ms=11)
    ax.text(4.52, 3.28, "否", fontsize=10, color=GRAY)

    ax.add_patch(FancyBboxPatch((0.35, 0.25), 11.3, 1.20,
                                boxstyle="round,pad=0.05,rounding_size=0.10",
                                fc="#FFF3E0", ec="#EF6C00", lw=1.3))
    ax.text(6.0, 1.25,
            "能降维：CPS（平面应力）/ CPE（平面应变）/ CAX（轴对称）",
            ha="center", va="center", fontsize=9.6, color="#22303A")
    ax.text(6.0, 0.90,
            "物理扩展：C3D8P、CPE4P（孔压）/ COH3D8、COH2D4（裂缝）",
            ha="center", va="center", fontsize=9.6, color="#22303A")
    ax.text(6.0, 0.55,
            "XFEM（路径未知）/ CEL（大变形、盐构造）/ B31 或 PIPE（管线）",
            ha="center", va="center", fontsize=9.6, color="#22303A")

    fig.text(0.5, 0.975, "石油工程常见选型：从几何到单元名",
             fontsize=15.0, ha="center", va="center")
    save(fig, "fig40-element-selection.png")


if __name__ == "__main__":
    fig01()
    fig02()
    fig03()
    fig04()
    fig05()
    fig06()
    fig07()
    fig08()
    fig09()
    fig10()
    fig11()
    fig12()
    fig13()
    fig14()
    fig15()
    fig16()
    fig17()
    fig18()
    fig19()
    fig20()
    fig21()
    fig22()
    fig23()
    fig24()
    fig25()
    fig26()
    fig27()
    fig28()
    fig29()
    fig30()
    fig31()
    fig32()
    fig33()
    fig34()
    fig35()
    fig36()
    fig37()
    fig38()
    fig39()
    fig40()
