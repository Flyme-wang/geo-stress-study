#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""单元类型笔记的参考实现（二维平面应力）。

本模块不是通用有限元程序，只服务于 docs/06-numerics/element-types.md 的
可核对数值：

- Q4（双线性四边形）、Q8（8 节点 Serendipity）、T3（常应变三角形）
  的形函数、B 矩阵与刚度矩阵；
- 2x2 / 1x1 / 3x3 高斯积分；
- 悬臂梁网格、端部集中荷载与挠度求解；
- 单元特征值谱（刚体模态、沙漏模态计数）；
- 纯弯曲下的寄生剪切（剪切锁死）解析值与数值值；
- 辐条屈曲、网格/自由度账、平面应力与平面应变的 sigma_z 换算。

单位约定：长度 mm，模量 MPa，力 N，应力 MPa。
所有函数只依赖 numpy。
"""

import numpy as np

# ---------------------------------------------------------------- 基准算例
E_STEEL = 210000.0      # MPa
NU_STEEL = 0.3
THICKNESS = 10.0        # mm，平面应力厚度
BEAM_L = 100.0          # mm
BEAM_H = 10.0           # mm
TIP_LOAD = 1000.0       # N


def plane_stress_D(E=E_STEEL, nu=NU_STEEL):
    """各向同性平面应力弹性矩阵（3x3）。"""
    return E / (1.0 - nu ** 2) * np.array([
        [1.0, nu, 0.0],
        [nu, 1.0, 0.0],
        [0.0, 0.0, (1.0 - nu) / 2.0],
    ])


# ---------------------------------------------------------------- 形函数
def q4_shape(xi, eta):
    """双线性四边形：节点顺序 (0,0) (1,0) (1,1) (0,1)。

    返回 N（4）与 dN/dxi、dN/deta（4x2）。
    """
    N = 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ])
    dN = 0.25 * np.array([
        [-(1 - eta), -(1 - xi)],
        [(1 - eta), -(1 + xi)],
        [(1 + eta), (1 + xi)],
        [-(1 + eta), (1 - xi)],
    ])
    return N, dN


def q8_shape(xi, eta):
    """8 节点 Serendipity 四边形：4 角点 + 4 边中点。

    节点顺序：左下、右下、右上、左上、下中、右中、上中、左中。
    """
    N = np.array([
        -0.25 * (1 - xi) * (1 - eta) * (1 + xi + eta),
        -0.25 * (1 + xi) * (1 - eta) * (1 - xi + eta),
        -0.25 * (1 + xi) * (1 + eta) * (1 - xi - eta),
        -0.25 * (1 - xi) * (1 + eta) * (1 + xi - eta),
        0.5 * (1 - xi ** 2) * (1 - eta),
        0.5 * (1 + xi) * (1 - eta ** 2),
        0.5 * (1 - xi ** 2) * (1 + eta),
        0.5 * (1 - xi) * (1 - eta ** 2),
    ])
    dN = np.array([
        [0.25 * (1 - eta) * (2 * xi + eta), 0.25 * (1 - xi) * (xi + 2 * eta)],
        [0.25 * (1 - eta) * (2 * xi - eta), 0.25 * (1 + xi) * (-xi + 2 * eta)],
        [0.25 * (1 + eta) * (2 * xi + eta), 0.25 * (1 + xi) * (xi + 2 * eta)],
        [0.25 * (1 + eta) * (2 * xi - eta), 0.25 * (1 - xi) * (-xi + 2 * eta)],
        [-xi * (1 - eta), -0.5 * (1 - xi ** 2)],
        [0.5 * (1 - eta ** 2), -(1 + xi) * eta],
        [-xi * (1 + eta), 0.5 * (1 - xi ** 2)],
        [-0.5 * (1 - eta ** 2), -(1 - xi) * eta],
    ])
    return N, dN


def q8_unit_square():
    """单位正方形的 Q8 节点坐标，与 q8_shape 的顺序一致。"""
    return np.array([
        [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0],
        [0.5, 0.0], [1.0, 0.5], [0.5, 1.0], [0.0, 0.5],
    ])


# ---------------------------------------------------------------- 积分
def gauss_2d(n):
    """二维张量积高斯积分：返回 [(xi, eta, w), ...]。

    n=1 为减缩积分，n=2 为常用全积分，n=3 用于二次单元。
    """
    if n == 1:
        x, w = [0.0], [2.0]
    elif n == 2:
        g = 1.0 / np.sqrt(3.0)
        x, w = [-g, g], [1.0, 1.0]
    elif n == 3:
        g = np.sqrt(3.0 / 5.0)
        x, w = [-g, 0.0, g], [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
    else:
        raise ValueError("只实现 1/2/3 点高斯积分")
    return [(x[i], x[j], w[i] * w[j])
            for i in range(len(x)) for j in range(len(x))]


# ---------------------------------------------------------------- 刚度
def isoparametric_stiffness(coords, shape, gauss, D=None, thickness=THICKNESS):
    """等参单元平面应力刚度：K = sum B^T D B detJ w t。"""
    if D is None:
        D = plane_stress_D()
    coords = np.asarray(coords, float)
    n = len(coords)
    K = np.zeros((2 * n, 2 * n))
    for xi, eta, w in gauss:
        _, dN = shape(xi, eta)
        J = dN.T @ coords
        detJ = np.linalg.det(J)
        if detJ <= 0.0:
            raise ValueError("detJ <= 0：节点顺序或坐标有误")
        dNdx = dN @ np.linalg.inv(J)
        B = np.zeros((3, 2 * n))
        B[0, 0::2] = dNdx[:, 0]
        B[1, 1::2] = dNdx[:, 1]
        B[2, 0::2] = dNdx[:, 1]
        B[2, 1::2] = dNdx[:, 0]
        K += (B.T @ D @ B) * detJ * w * thickness
    return K


def cst_stiffness(coords, D=None, thickness=THICKNESS):
    """常应变三角形（T3/CST）刚度；三节点逆时针。"""
    if D is None:
        D = plane_stress_D()
    coords = np.asarray(coords, float)
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]
    area2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if area2 <= 0.0:
        raise ValueError("三角形面积为非正：检查节点顺序")
    B = np.array([
        [y2 - y3, 0.0, y3 - y1, 0.0, y1 - y2, 0.0],
        [0.0, x3 - x2, 0.0, x1 - x3, 0.0, x2 - x1],
        [x3 - x2, y2 - y3, x1 - x3, y3 - y1, x2 - x1, y1 - y2],
    ]) / area2
    return (B.T @ D @ B) * (0.5 * area2) * thickness


# ---------------------------------------------------------------- 网格
def mesh_q4(nx, ny, L=BEAM_L, H=BEAM_H):
    """规则 Q4 网格：nx 段沿 x，ny 段沿 y。"""
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(0.0, H, ny + 1)
    nodes = np.array([(x, y) for y in ys for x in xs])
    idx = lambda i, j: j * (nx + 1) + i
    elems = []
    for j in range(ny):
        for i in range(nx):
            elems.append([idx(i, j), idx(i + 1, j),
                          idx(i + 1, j + 1), idx(i, j + 1)])
    return nodes, elems


def mesh_q8(nx, ny, L=BEAM_L, H=BEAM_H):
    """规则 Q8 网格：每条 Q4 边中点补一个节点。"""
    xs = np.linspace(0.0, L, 2 * nx + 1)
    ys = np.linspace(0.0, H, 2 * ny + 1)
    all_nodes = [(x, y) for y in ys for x in xs]
    pos = {p: k for k, p in enumerate(all_nodes)}
    cnodes, quads = mesh_q4(nx, ny, L, H)
    elems = []
    for e in quads:
        corners = [cnodes[j] for j in e]
        ij = [(int(round(x / L * 2 * nx)), int(round(y / H * 2 * ny)))
              for x, y in corners]
        i0, j0 = ij[0]
        mids = [(i0 + 1, j0), (i0 + 2, j0 + 1),
                (i0 + 1, j0 + 2), (i0, j0 + 1)]
        nn = [pos[(xs[i0], ys[j0])], pos[(xs[i0 + 2], ys[j0])],
              pos[(xs[i0 + 2], ys[j0 + 2])], pos[(xs[i0], ys[j0 + 2])]]
        nn += [pos[(xs[i], ys[j])] for i, j in mids]
        elems.append(nn)
    used = sorted({j for e in elems for j in e})
    remap = {j: k for k, j in enumerate(used)}
    nodes = np.array([all_nodes[j] for j in used])
    elems = [[remap[j] for j in e] for e in elems]
    return nodes, elems


def mesh_t3(nx, ny, L=BEAM_L, H=BEAM_H):
    """每个 Q4 沿对角线切成两个 T3。"""
    nodes, quads = mesh_q4(nx, ny, L, H)
    elems = []
    for a, b, c, d in quads:
        elems.append([a, b, c])
        elems.append([a, c, d])
    return nodes, elems


# ---------------------------------------------------------------- 求解
def assemble_and_solve(nodes, elems, stiffness_fn, load=TIP_LOAD,
                       L=BEAM_L):
    """左端全固定、右端竖向集中力的平面应力求解，返回全部节点位移。"""
    nodes = np.asarray(nodes, float)
    nn = len(nodes)
    K = np.zeros((2 * nn, 2 * nn))
    for e in elems:
        idx = []
        for i in e:
            idx += [2 * i, 2 * i + 1]
        K[np.ix_(idx, idx)] += stiffness_fn(nodes[e])

    F = np.zeros(2 * nn)
    tip = sorted([i for i, p in enumerate(nodes) if abs(p[0] - L) < 1e-9],
                 key=lambda i: nodes[i][1])
    shares = {2: [0.5, 0.5], 3: [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0]}
    w = shares.get(len(tip), [1.0 / len(tip)] * len(tip))
    for i, s in zip(tip, w):
        F[2 * i + 1] -= load * s

    fixed = []
    for i, p in enumerate(nodes):
        if abs(p[0]) < 1e-9:
            fixed += [2 * i, 2 * i + 1]
    free = np.array([i for i in range(2 * nn) if i not in fixed])
    U = np.zeros(2 * nn)
    U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    return U


def tip_deflection(nodes, elems, stiffness_fn, load=TIP_LOAD, L=BEAM_L):
    """端部平均竖向挠度，向下为正（mm）。"""
    nodes = np.asarray(nodes, float)
    U = assemble_and_solve(nodes, elems, stiffness_fn, load, L)
    tip = [i for i, p in enumerate(nodes) if abs(p[0] - L) < 1e-9]
    return -np.mean([U[2 * i + 1] for i in tip])


def beam_tip_exact(load=TIP_LOAD, L=BEAM_L, H=BEAM_H,
                   E=E_STEEL, thickness=THICKNESS):
    """Euler–Bernoulli 悬臂梁端部挠度 P L^3 / (3 E I)，向下为正。"""
    I = thickness * H ** 3 / 12.0
    return load * L ** 3 / (3.0 * E * I)


def centerline(nodes, U, L=BEAM_L, tol=1e-9):
    """按 x 坐标聚合节点，返回 (x, 平均竖向位移)；用于画挠曲线。"""
    nodes = np.asarray(nodes, float)
    xs = sorted({round(float(p[0]), 9) for p in nodes
                 if abs(p[0] - L) < L + 1.0})
    out_x, out_u = [], []
    for x in xs:
        ids = [i for i, p in enumerate(nodes) if abs(p[0] - x) < tol]
        out_x.append(x)
        out_u.append(np.mean([U[2 * i + 1] for i in ids]))
    return np.array(out_x), -np.array(out_u)


def beam_centerline_exact(load=TIP_LOAD, L=BEAM_L, H=BEAM_H,
                          E=E_STEEL, thickness=THICKNESS, n=200):
    """Euler–Bernoulli 挠曲线 delta(x) = P x^2 (3L - x) / (6 E I)。"""
    I = thickness * H ** 3 / 12.0
    x = np.linspace(0.0, L, n)
    return x, load * x ** 2 * (3.0 * L - x) / (6.0 * E * I)


# ---------------------------------------------------------------- 特征值
def element_eigenvalues(shape, gauss, D=None, thickness=1.0):
    """单位正方形单元的刚度特征值，升序返回。"""
    if shape is q8_shape:
        coords = q8_unit_square()
    else:
        coords = np.array([[0.0, 0.0], [1.0, 0.0],
                           [1.0, 1.0], [0.0, 1.0]])
    K = isoparametric_stiffness(coords, shape, gauss, D, thickness)
    return np.sort(np.linalg.eigvalsh(K))


def zero_mode_count(eigs, rel=1e-8):
    """小于最大特征值 rel 倍的特征值个数（刚体 + 沙漏）。"""
    return int(np.sum(np.abs(eigs) < rel * np.max(np.abs(eigs))))


def hourglass_modes():
    """Q4 减缩积分单刚的 2 个沙漏模态（已剔除 3 个刚体模态）。

    返回 2x8 数组，每行一个模态，按最大分量归一化。
    """
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    K = isoparametric_stiffness(coords, q4_shape, gauss_2d(1))
    w, V = np.linalg.eigh(K)
    null = V[:, np.argsort(np.abs(w))[:5]]
    R = np.zeros((8, 3))
    for k, (x, y) in enumerate(coords):
        R[2 * k, 0] = 1.0
        R[2 * k + 1, 1] = 1.0
        R[2 * k, 2] = -y
        R[2 * k + 1, 2] = x
    P = null - R @ np.linalg.lstsq(R, null, rcond=None)[0]
    u, _, _ = np.linalg.svd(P)
    modes = u[:, :2].T
    for m in modes:
        m /= np.max(np.abs(m))
    return modes


# ---------------------------------------------------------------- 剪切锁死
def parasitic_shear(kappa=1.0, a=10.0):
    """纯弯曲下 Q4 的寄生剪应变 gamma = kappa a / (2 sqrt(3))。

    与单元厚度 h 无关；与弯曲正应变 kappa h 相比，量级为 a/h。
    """
    return kappa * a / (2.0 * np.sqrt(3.0))


def pure_bending_strains(a=10.0, h=1.0, kappa=1.0):
    """在 2x2 高斯点上量 Q4 的 strain：返回 [(eps_x, gamma_xy), ...]。

    节点位移取自 u_x = -kappa x y，u_y = kappa x^2 / 2，y 以中面为 0，
    即单元占据 x∈[0,a]、y∈[-h/2,h/2]，中性轴在 y=0。
    """
    coords = np.array([[0.0, -h / 2], [a, -h / 2], [a, h / 2], [0.0, h / 2]])
    ue = np.zeros(8)
    for k, (x, y) in enumerate(coords):
        ue[2 * k] = -kappa * x * y
        ue[2 * k + 1] = 0.5 * kappa * x ** 2
    out = []
    for xi, eta, _ in gauss_2d(2):
        _, dN = q4_shape(xi, eta)
        J = dN.T @ coords
        dNdx = dN @ np.linalg.inv(J)
        B = np.zeros((3, 8))
        B[0, 0::2] = dNdx[:, 0]
        B[1, 1::2] = dNdx[:, 1]
        B[2, 0::2] = dNdx[:, 1]
        B[2, 1::2] = dNdx[:, 0]
        e = B @ ue
        out.append((e[0], e[2]))
    return out


def q8_pure_bending_shear(a=10.0, h=1.0, kappa=1.0):
    """Q8 在同样的纯弯曲位移场下，所有 3x3 高斯点的剪应变最大值。"""
    coords = np.array([[0.0, 0.0], [a, 0.0], [a, h], [0.0, h],
                       [a / 2, 0.0], [a, h / 2], [a / 2, h], [0.0, h / 2]])
    ue = np.zeros(16)
    for k, (x, y) in enumerate(coords):
        ue[2 * k] = -kappa * x * y
        ue[2 * k + 1] = 0.5 * kappa * x ** 2
    worst = 0.0
    for xi, eta, _ in gauss_2d(3):
        _, dN = q8_shape(xi, eta)
        J = dN.T @ coords
        dNdx = dN @ np.linalg.inv(J)
        B = np.zeros((3, 16))
        B[0, 0::2] = dNdx[:, 0]
        B[1, 1::2] = dNdx[:, 1]
        B[2, 0::2] = dNdx[:, 1]
        B[2, 1::2] = dNdx[:, 0]
        worst = max(worst, abs((B @ ue)[2]))
    return worst


# ---------------------------------------------------------------- 工程账
def dof_accounting():
    """原文挡泥板算例：实体 1 mm 立方 vs 壳 5 mm 中面网格。"""
    solid = {"nx": 600, "ny": 100, "nz": 4}
    solid["elements"] = solid["nx"] * solid["ny"] * solid["nz"]
    solid["nodes"] = ((solid["nx"] + 1) * (solid["ny"] + 1)
                      * (solid["nz"] + 1))
    solid["dof"] = 3 * solid["nodes"]

    shell = {"nx": 600 // 5, "ny": 100 // 5}
    shell["elements"] = shell["nx"] * shell["ny"]
    shell["nodes"] = (shell["nx"] + 1) * (shell["ny"] + 1)
    shell["dof"] = 6 * shell["nodes"]
    return {"solid": solid, "shell": shell}


def spoke_buckling(d=2.0, L=290.0, E=E_STEEL):
    """辐条 Euler 屈曲荷载（两端铰支）：Pcr = pi^2 E I / L^2。"""
    I = np.pi * d ** 4 / 64.0
    return np.pi ** 2 * E * I / L ** 2


def spoke_tension_stress(mass_kg=100.0, d=2.0, g=9.80665):
    """辐条拉 100 kg 时的名义拉应力（MPa）。"""
    area = np.pi * (d / 2.0) ** 2
    return mass_kg * g / area


def plane_strain_sigma_z(sx, sy, nu=NU_STEEL):
    """平面应变 epsilon_z = 0 时的 sigma_z = nu (sigma_x + sigma_y)。"""
    return nu * (sx + sy)


def plane_stress_eps_z(sx, sy, E=E_STEEL, nu=NU_STEEL):
    """平面应力 sigma_z = 0 时的 epsilon_z = -nu (sigma_x + sigma_y) / E。"""
    return -nu * (sx + sy) / E
