"""CPS3 完全拉格朗日（Total Lagrangian, TL）教学参考实现。

符号与 M6-10 笔记一致：

- 参考坐标 $(X, Y)$，位移 $(u, v)$；Voigt 应变
  $E = [E_{XX}, E_{YY}, 2E_{XY}]^T$，剪切分量是**工程剪应变**；
- 第二 Piola–Kirchhoff 应力 $S = [S_{XX}, S_{YY}, S_{XY}]^T$，
  剪切分量是**张量剪应力**，与工程剪应变点积后正好等于张量双点积；
- 单元内力 $f_{int} = V_0 B^T S$，切线 $K = K_{mat} + K_{geo}$；
- 所有积分在初始参考构形上完成（TL），位移梯度与应变为当前构形量。

这个模块只服务教学与数值校验，不追求性能，也不处理接触、塑性或孔压。
"""

from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------
# 几何
# --------------------------------------------------------------------------
def tri_geometry(coords, thickness=1.0):
    """返回参考构形中三角形单元的几何量。

    Parameters
    ----------
    coords : (3, 2) array
        三个节点的初始坐标，必须逆时针编号，$A_0 > 0$。
    thickness : float
        均匀厚度 $h_0$。

    Returns
    -------
    dict with keys ``A``, ``V``, ``h``, ``b``, ``c``, ``dNdx``, ``dNdy``.
    """
    coords = np.asarray(coords, dtype=float)
    if coords.shape != (3, 2):
        raise ValueError("coords 必须是 3x2 数组")
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]

    two_a = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if two_a <= 0.0:
        raise ValueError("节点须逆时针编号，使 A0 > 0")
    area = 0.5 * two_a

    # N_a = (a_a + b_a X + c_a Y) / (2 A0)
    b = np.array([y2 - y3, y3 - y1, y1 - y2], dtype=float)
    c = np.array([x3 - x2, x1 - x3, x2 - x1], dtype=float)

    return {
        "coords": coords,
        "A": area,
        "V": area * thickness,
        "h": float(thickness),
        "b": b,
        "c": c,
        "dNdx": b / (2.0 * area),
        "dNdy": c / (2.0 * area),
    }


def shape_functions(geom, points):
    """在参考坐标点上计算三个线性形函数，返回 (n_points, 3)。"""
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    a = np.array([
        geom["coords"][1, 0] * geom["coords"][2, 1]
        - geom["coords"][2, 0] * geom["coords"][1, 1],
        geom["coords"][2, 0] * geom["coords"][0, 1]
        - geom["coords"][0, 0] * geom["coords"][2, 1],
        geom["coords"][0, 0] * geom["coords"][1, 1]
        - geom["coords"][1, 0] * geom["coords"][0, 1],
    ], dtype=float)
    return (a[None, :]
            + pts[:, 0:1] * geom["b"][None, :]
            + pts[:, 1:2] * geom["c"][None, :]) / (2.0 * geom["A"])


def gradient_matrix(geom):
    """位移梯度映射 $\\theta = G d$，$G$ 为 4x6，只由参考几何决定。"""
    g = np.zeros((4, 6))
    for a in range(3):
        g[0, 2 * a] = geom["dNdx"][a]
        g[1, 2 * a] = geom["dNdy"][a]
        g[2, 2 * a + 1] = geom["dNdx"][a]
        g[3, 2 * a + 1] = geom["dNdy"][a]
    return g


# --------------------------------------------------------------------------
# 运动学
# --------------------------------------------------------------------------
def displacement_gradient(geom, d):
    """返回 $\\theta = [u_{,X}, u_{,Y}, v_{,X}, v_{,Y}]^T$ 与 $G$。"""
    return gradient_matrix(geom) @ np.asarray(d, dtype=float).reshape(6)


def deformation_gradient(geom, d):
    """返回 $F = I + \\partial u / \\partial X$（2x2）。"""
    theta = displacement_gradient(geom, d)
    return np.array([[1.0 + theta[0], theta[1]],
                     [theta[2], 1.0 + theta[3]]])


def gl_strain(geom, d):
    """Green–Lagrange 应变（Voigt，工程剪应变）：$[E_{XX}, E_{YY}, 2E_{XY}]^T$。"""
    t1, t2, t3, t4 = displacement_gradient(geom, d)
    return np.array([
        t1 + 0.5 * (t1 * t1 + t3 * t3),
        t4 + 0.5 * (t2 * t2 + t4 * t4),
        t2 + t3 + t1 * t2 + t3 * t4,
    ])


def small_strain(geom, d):
    """小应变（Voigt，工程剪应变），用于与 GL 应变对照。"""
    t1, t2, t3, t4 = displacement_gradient(geom, d)
    return np.array([t1, t4, t2 + t3])


def b_matrix(geom, d):
    """全几何 B 矩阵（3x6）：$\\delta E = B(d)\\,\\delta d$，B = B_L + B_NL。"""
    t1, t2, t3, t4 = displacement_gradient(geom, d)
    b = np.zeros((3, 6))
    for a in range(3):
        nx = geom["dNdx"][a]
        ny = geom["dNdy"][a]
        # 线性部分（与小应变 CST 相同）
        b[0, 2 * a] = nx
        b[1, 2 * a + 1] = ny
        b[2, 2 * a] = ny
        b[2, 2 * a + 1] = nx
        # 非线性部分（随当前位移梯度变化）
        b[0, 2 * a] += t1 * nx
        b[0, 2 * a + 1] += t3 * nx
        b[1, 2 * a] += t2 * ny
        b[1, 2 * a + 1] += t4 * ny
        b[2, 2 * a] += t1 * ny + t2 * nx
        b[2, 2 * a + 1] += t3 * ny + t4 * nx
    return b


def b_linear(geom):
    """线性 B 矩阵 B_L（CST 应变矩阵，3x6）。"""
    return b_matrix(geom, np.zeros(6))


def b_nonlinear(geom, d):
    """非线性 B 矩阵 B_NL(d) = B(d) - B_L。"""
    return b_matrix(geom, d) - b_linear(geom)


# --------------------------------------------------------------------------
# 本构：平面应力 St. Venant–Kirchhoff
# --------------------------------------------------------------------------
def stvk_plane_stress_moduli(young, poisson):
    """平面应力 StVK 的 Voigt 材料矩阵 C（3x3）。

    应变向量含工程剪应变 $2E_{XY}$，应力向量含张量剪应力 $S_{XY}$，
    所以剪切对角项是 $(1-\\nu)/2$，不是 $(1-2\\nu)/2$。
    """
    e = float(young)
    nu = float(poisson)
    if not (-1.0 < nu < 0.5):
        raise ValueError("泊松比须在 (-1, 0.5) 内")
    return e / (1.0 - nu * nu) * np.array([
        [1.0, nu, 0.0],
        [nu, 1.0, 0.0],
        [0.0, 0.0, 0.5 * (1.0 - nu)],
    ])


def stvk_stress(strain_voigt, young, poisson):
    """第二 PK 应力 $S = C E$（Voigt）。"""
    return stvk_plane_stress_moduli(young, poisson) @ np.asarray(
        strain_voigt, dtype=float)


# --------------------------------------------------------------------------
# 内力与切线
# --------------------------------------------------------------------------
def internal_force(geom, d, young, poisson):
    """单元内力 $f_{int} = V_0 B^T S$（6，向量）。"""
    strain = gl_strain(geom, d)
    stress = stvk_stress(strain, young, poisson)
    return geom["V"] * b_matrix(geom, d).T @ stress, strain, stress


def tangent(geom, d, young, poisson):
    """返回 (K_mat, K_geo, K_total)，均为 6x6。"""
    d = np.asarray(d, dtype=float).reshape(6)
    strain = gl_strain(geom, d)
    stress = stvk_stress(strain, young, poisson)
    sxx, syy, sxy = stress
    b = b_matrix(geom, d)
    c = stvk_plane_stress_moduli(young, poisson)

    k_mat = geom["V"] * b.T @ c @ b

    # K_geo = V0 G^T diag(S0, S0) G，S0 = [[Sxx, Sxy], [Sxy, Syy]]
    g = gradient_matrix(geom)
    s0 = np.array([[sxx, sxy], [sxy, syy]])
    m = np.zeros((4, 4))
    m[0:2, 0:2] = s0
    m[2:4, 2:4] = s0
    k_geo = geom["V"] * g.T @ m @ g

    return k_mat, k_geo, k_mat + k_geo


# --------------------------------------------------------------------------
# 两单元 patch 与 Newton
# --------------------------------------------------------------------------
def square_patch():
    """单位正方形、两个 CPS3 单元（逆时针）。

    节点编号：0 (0,0)、1 (1,0)、2 (1,1)、3 (0,1)；
    单元：T1 = [0, 1, 2]，T2 = [0, 2, 3]。
    """
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    conn = [(0, 1, 2), (0, 2, 3)]
    return coords, conn


def assemble_patch(coords, conn, d, young, poisson, thickness=1.0,
                   material_only=False):
    """组装 patch 的内力（8）与切线（8x8）。"""
    d = np.asarray(d, dtype=float).reshape(8)
    f_int = np.zeros(8)
    k = np.zeros((8, 8))
    for tri in conn:
        geom = tri_geometry(coords[list(tri)], thickness)
        de = d[[2 * n + i for n in tri for i in range(2)]]
        fe, _, _ = internal_force(geom, de, young, poisson)
        k_mat, k_geo, k_tot = tangent(geom, de, young, poisson)
        ke = k_mat if material_only else k_tot
        idx = [2 * n + i for n in tri for i in range(2)]
        for i_local, i_global in enumerate(idx):
            f_int[i_global] += fe[i_local]
            for j_local, j_global in enumerate(idx):
                k[i_global, j_global] += ke[i_local, j_local]
    return f_int, k


def newton_patch(load, young, poisson, material_only=False, max_iter=30,
                 tol=1e-12, thickness=1.0):
    """左边界固定、右边界受均布死载荷的 patch Newton 求解。

    Parameters
    ----------
    load : float
        右边界单位参考长度的总拉力（厚度 1 时即合力）。

    Returns
    -------
    dict: d, iterations, residuals, converged, f_ext
    """
    coords, conn = square_patch()
    d = np.zeros(8)
    # 自由度：0,1 = 节点 0；2,3 = 节点 1；4,5 = 节点 2；6,7 = 节点 3
    fixed = [0, 1, 6, 7]
    free = [i for i in range(8) if i not in fixed]
    f_ext = np.zeros(8)
    f_ext[2] += 0.5 * load
    f_ext[4] += 0.5 * load

    residuals = []
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        f_int, k = assemble_patch(coords, conn, d, young, poisson,
                                  thickness, material_only)
        r = f_ext - f_int
        r_free = r[free]
        norm = np.linalg.norm(r_free)
        residuals.append(norm)
        if norm <= tol * max(1.0, np.linalg.norm(f_ext)):
            converged = True
            break
        k_ff = k[np.ix_(free, free)]
        delta = np.linalg.solve(k_ff, r_free)
        d[free] += delta
    return {
        "d": d,
        "iterations": it,
        "residuals": residuals,
        "converged": converged,
        "f_ext": f_ext,
        "coords": coords,
        "conn": conn,
    }


def uniform_stretch_state(stretch, young, poisson, thickness=1.0):
    """右手直角三角形的均匀拉伸状态，供解析核对。

    三角形 (0,0)、(2,0)、(0,1)，位移 $u = sX, v = 0$。
    """
    coords = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 1.0]])
    geom = tri_geometry(coords, thickness)
    d = np.array([stretch * 0.0, 0.0,
                  stretch * 2.0, 0.0,
                  stretch * 0.0, 0.0])
    strain = gl_strain(geom, d)
    stress = stvk_stress(strain, young, poisson)
    f, _, _ = internal_force(geom, d, young, poisson)
    return geom, d, strain, stress, f