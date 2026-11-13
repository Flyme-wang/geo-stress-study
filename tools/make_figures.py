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
                                FancyBboxPatch)

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
