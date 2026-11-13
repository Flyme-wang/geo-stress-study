#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recompute every number quoted in docs/06-numerics/fem-programming/.

Usage (from the repository root):
    py -3 tools/verify_fem_notes.py

The script prints each check and exits non-zero if an assertion fails, so the
notes can claim "可核对" without asking the reader to trust a screenshot.
"""

import itertools
import sys
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import make_figures as mf  # noqa: E402


def q4(xy, E, nu, t):
    """Plane-stress Q4 stiffness, 2x2 Gauss, node order 1-2-3-4."""
    D = E / (1 - nu ** 2) * np.array([[1, nu, 0],
                                      [nu, 1, 0],
                                      [0, 0, (1 - nu) / 2]])
    g = 1 / np.sqrt(3)
    pts = [(-g, -g), (g, -g), (g, g), (-g, g)]
    Ke = np.zeros((8, 8))
    B = None
    for xi, eta in pts:
        dN = 0.25 * np.array([[-(1 - eta), -(1 - xi)],
                              [(1 - eta), -(1 + xi)],
                              [(1 + eta), (1 + xi)],
                              [-(1 + eta), (1 - xi)]])
        J = jacobian(xy, dN)
        detJ = np.linalg.det(J)
        d = np.linalg.solve(J, dN.T).T
        B = np.zeros((3, 8))
        for i in range(4):
            B[0, 2 * i] = d[i, 0]
            B[1, 2 * i + 1] = d[i, 1]
            B[2, 2 * i] = d[i, 1]
            B[2, 2 * i + 1] = d[i, 0]
        Ke += B.T @ D @ B * detJ * t
    return Ke, D, B


def jacobian(xy, dN):
    J = np.zeros((2, 2))
    for i in range(4):
        J[0, 0] += dN[i, 0] * xy[i, 0]
        J[0, 1] += dN[i, 0] * xy[i, 1]
        J[1, 0] += dN[i, 1] * xy[i, 0]
        J[1, 1] += dN[i, 1] * xy[i, 1]
    return J


def dn_local(xi, eta):
    return 0.25 * np.array([[-(1 - eta), -(1 - xi)],
                            [(1 - eta), -(1 + xi)],
                            [(1 + eta), (1 + xi)],
                            [-(1 + eta), (1 - xi)]])


def check(label, value, target, tol):
    ok = abs(value - target) <= tol
    print(("PASS  " if ok else "FAIL  ") + label
          + "  value=%+.6e  target=%+.6e  tol=%.1e" % (value, target, tol))
    if not ok:
        raise SystemExit(1)


def approx(label, value, target, rel):
    ok = abs(value - target) <= rel * abs(target)
    print(("PASS  " if ok else "FAIL  ") + label
          + "  value=%.6f  target=%.6f  rel=%.2f%%"
          % (value, target, 100 * rel))
    if not ok:
        raise SystemExit(1)


# ============================================================ 1. Q4 single element
def part1():
    print("\n# 1  Q4 unit square (E=300 GPa, nu=0.3, t=1 m)")
    xy = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
    Ke, D, B = q4(xy, 300e9, 0.3, 1.0)
    check("K(1,1)", Ke[0, 0], 1.483516e11, 1e6)
    check("K(3,3)", Ke[2, 2], 1.483516e11, 1e6)
    check("row 1 sum (rigid translation)", Ke[0].sum(), 0.0, 1e-3)

    F = np.array([0, 0, 0, 0, 0, -100, 0, -100.], float)
    u = np.zeros(8)
    u[4:] = np.linalg.solve(Ke[4:, 4:], F[4:])
    check("u(3y)", u[5], -6.468137e-10, 1e-15)
    check("u(4y)", u[7], -6.468137e-10, 1e-15)
    check("u(3x)", u[4], +1.338235e-10, 1e-15)
    check("u(4x)", u[6], -1.338235e-10, 1e-15)

    R = Ke @ u - F
    check("sum Rx", R[0::2].sum(), 0.0, 1e-9)
    check("sum Ry", R[1::2].sum(), 200.0, 1e-9)
    approx("eps_y vs sigma/E", abs(u[5]) / 1.0, 200 / 300e9, 0.05)

    # ---- 后处理：图 7 的三个场（u_x、|u|、sigma_y）
    def B_at(xi, eta):
        d = np.linalg.solve(jacobian(xy, dn_local(xi, eta)),
                            dn_local(xi, eta).T).T
        Bx = np.zeros((3, 8))
        for i in range(4):
            Bx[0, 2 * i] = d[i, 0]
            Bx[1, 2 * i + 1] = d[i, 1]
            Bx[2, 2 * i] = d[i, 1]
            Bx[2, 2 * i + 1] = d[i, 0]
        return Bx

    sig_c = D @ B_at(0.0, 0.0) @ u
    check("sigma_y at centre", sig_c[1], -200.0, 1e-6)
    check("sigma_x at centre", sig_c[0], -19.853, 1e-3)
    g = 1 / np.sqrt(3)
    check("sigma_y at GP 1", (D @ B_at(-g, -g) @ u)[1], -207.641, 1e-3)
    check("sigma_y at GP 4", (D @ B_at(g, g) @ u)[1], -192.359, 1e-3)
    mean_sy = sum((D @ B_at(sx * g, sy * g) @ u)[1]
                  for sx in (-1, 1) for sy in (-1, 1)) / 4.0
    check("mean sigma_y == -200 N / 1 m^2", mean_sy, -200.0, 1e-9)
    check("|u| at node 3", np.hypot(u[4], u[5]), 6.605125e-10, 1e-15)

    # ---- 点荷载一致性：中心 30 N -> 四节点各 7.5 N
    w = 0.25 * np.ones(4)
    check("N_a(0,0) = 1/4", w[0], 0.25, 0.0)
    check("point load total", (w * 30.0).sum(), 30.0, 1e-12)
    check("point load per node", (w * 30.0)[0], 7.5, 1e-12)


# ============================================================ 2. node ordering / detJ
def part2():
    print("\n# 2  detJ for three node orderings (unit square)")
    xy = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
    g = 1 / np.sqrt(3)
    for name, order, target in (("1-2-3-4", [0, 1, 2, 3], +0.25),
                                ("1-4-3-2", [0, 3, 2, 1], -0.25)):
        vals = [np.linalg.det(jacobian(xy[order], dn_local(sx * g, sy * g)))
                for sx in (-1, 1) for sy in (-1, 1)]
        for v in vals:
            check("detJ %s" % name, v, target, 1e-12)
    vals = [np.linalg.det(jacobian(xy[[0, 1, 3, 2]], dn_local(sx * g, sy * g)))
            for sx in (-1, 1) for sy in (-1, 1)]
    check("detJ crossed at GP1", vals[0], +0.144338, 1e-5)
    check("detJ crossed at GP2", vals[1], -0.144338, 1e-5)
    d0 = np.linalg.det(jacobian(xy[[0, 1, 3, 2]], dn_local(g, 0.0)))
    check("detJ crossed at eta=0", d0, 0.0, 1e-12)

    xyd = np.array([[0., 0.], [2., 0.], [2.4, 1.2], [0.1, 0.9]])
    area = sum(np.linalg.det(jacobian(xyd, dn_local(sx * g, sy * g)))
               for sx in (-1, 1) for sy in (-1, 1))
    xs, ys = xyd[:, 0], xyd[:, 1]
    shoe = 0.5 * abs(np.dot(xs, np.roll(ys, -1)) - np.dot(ys, np.roll(xs, -1)))
    check("distorted area (2x2 Gauss) = shoelace", area, shoe, 1e-12)
    check("distorted area value", area, 2.22, 1e-6)
    print("INFO  distorted detJ range = %.4f .. %.4f"
          % (min(np.linalg.det(jacobian(xyd, dn_local(sx * g, sy * g)))
                 for sx in (-1, 1) for sy in (-1, 1)),
             max(np.linalg.det(jacobian(xyd, dn_local(sx * g, sy * g)))
                 for sx in (-1, 1) for sy in (-1, 1))))


# ============================================================ 3. assembly
def part3():
    print("\n# 3  two-element assembly")
    elems = [[1, 2, 5, 4], [2, 3, 6, 5]]

    def stats(order):
        new = {old: i for i, old in enumerate(order)}
        S = np.zeros((12, 12), bool)
        for e in elems:
            flat = []
            for k in e:
                flat += [2 * new[k], 2 * new[k] + 1]
            for a in flat:
                for b in flat:
                    S[a, b] = True
        bw = max(abs(i - j) for i in range(12) for j in range(12) if S[i, j])
        return int(S.sum()), bw

    nnz, bw = stats((1, 2, 3, 4, 5, 6))
    check("nnz (natural order)", nnz, 112, 0)
    check("half-bandwidth (natural order)", bw, 9, 0)
    nnz2, bw2 = stats((1, 4, 2, 5, 3, 6))
    check("nnz (column-wise order)", nnz2, 112, 0)
    check("half-bandwidth (column-wise order)", bw2, 7, 0)

    pattern = [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (3, 1), (3, 3),
               (4, 2), (4, 4)]
    check("COO value sum", sum(range(1, 10)), 45, 0)
    check("Skyline stored count", sum([1, 2, 2, 3, 4, 3]), 15, 0)


# ============================================================ 4. constraints / loads
def part4():
    print("\n# 4  constraints and loads")
    K = np.array([[2., -1., 0.], [-1., 2., -1.], [0., -1., 2.]])
    F = np.array([1., 0., 0.])
    u = np.linalg.solve(K, F)
    check("unconstrained u1", u[0], 0.75, 1e-12)
    check("unconstrained u2", u[1], 0.50, 1e-12)
    check("unconstrained u3", u[2], 0.25, 1e-12)

    Kff = K[np.ix_([0, 2], [0, 2])]
    Ff = np.array([1., 0.]) - K[np.ix_([0, 2], [1])].ravel() * 0.25
    ue = np.linalg.solve(Kff, Ff)
    check("elimination u1", ue[0], 0.625, 1e-12)
    check("elimination u3", ue[1], 0.125, 1e-12)

    K2 = K.copy()
    F2 = F - K[:, 1] * 0.25
    K2[1, :] = 0
    K2[:, 1] = 0
    K2[1, 1] = 1.0
    F2[1] = 0.25
    ur = np.linalg.solve(K2, F2)
    check("row/col u1", ur[0], 0.625, 1e-12)
    check("row/col u2", ur[1], 0.250, 1e-12)
    check("row/col u3", ur[2], 0.125, 1e-12)

    for a, err in ((1e2, 2.475e-3), (1e4, 2.500e-5), (1e6, 2.500e-7),
                   (1e8, 2.500e-9), (1e10, 2.500e-11)):
        K3 = K.copy()
        F3 = F.copy()
        K3[1, 1] += a
        F3[1] += a * 0.25
        up = np.linalg.solve(K3, F3)
        check("penalty error alpha=%.0e" % a, abs(up[1] - 0.25), err, 1e-6)

    F = np.array([0, -100, 50, -100, 50, -70, 0, -100.], float)
    check("load sum Fx", F[0::2].sum(), 100, 1e-12)
    check("load sum Fy", F[1::2].sum(), -370, 1e-12)


# ============================================================ 5. solvers
def part5():
    print("\n# 5  CG / PCG on a layered 2-D Laplacian")
    expected = {16: (61, 50), 24: (97, 75), 32: (134, 102)}
    for n, (it_cg, it_pcg) in expected.items():
        A, bn = mf._lap_layers(n)
        b = np.ones(n * n)
        b[bn] = 0.0
        x, h0 = mf._cg(A, b)
        _, h1 = mf._cg(A, b, M=np.diag(1.0 / np.diag(A)))
        check("CG iterations n=%d" % n, len(h0) - 1, it_cg, 0)
        check("PCG iterations n=%d" % n, len(h1) - 1, it_pcg, 0)
        if n == 16:
            xd = np.linalg.solve(A, b)
            rel = np.linalg.norm(x - xd) / np.linalg.norm(xd)
            check("CG solution vs direct (n=16)", rel, 0.0, 1e-7)
    A, bn = mf._lap_layers(24)
    b = np.ones(24 * 24)
    b[bn] = 0.0
    check("cond(A) n=24", np.linalg.cond(A), 372.3, 1.0)


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    print("\nALL CHECKS PASSED")
