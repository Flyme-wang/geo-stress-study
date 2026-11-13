#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recompute every number quoted in docs/06-numerics/variational-approximation.md.

Usage (from the repository root):
    py -3 tools/verify_variational_note.py

The note works one model problem from beginning to end:

    -u'' = f  on (0,1),  u(0) = u(1) = 0,  u = sin(pi x),  f = pi^2 sin(pi x)

with the bilinear form A(u,v) = int u' v' and the load F(v) = int f v.
Every quantity quoted in the text is recomputed here: the Ritz/Galerkin linear
system, the energy identity A(u,u)-A(uN,uN)=||u-uN||_A^2 and the resulting lower
bound, the orthogonality of the error, best approximation, the O(h) / O(h^2)
rates, the nonsymmetric 2x2 counter-example, and the sparsity counts.
The script exits non-zero on the first mismatch.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss

_GP, _GW = leggauss(8)          # 8-point Gauss on [-1, 1]
A_EXACT = np.pi ** 2 / 2.0      # A(u,u) = int_0^1 (pi cos pi x)^2 dx


# ---------------------------------------------------------------- helpers
def check(label, value, target, tol):
    ok = abs(value - target) <= tol
    print(("PASS  " if ok else "FAIL  ") + label
          + "  value=%.10g  target=%.10g  tol=%.1e" % (value, target, tol))
    if not ok:
        raise SystemExit(1)


def check_true(label, cond, extra=""):
    print(("PASS  " if cond else "FAIL  ") + label + ("  " + extra if extra else ""))
    if not cond:
        raise SystemExit(1)


def quad(fun, nq=4000):
    """Composite 8-point Gauss on (0,1); fun takes an array, returns an array."""
    return quad_on(fun, np.linspace(0.0, 1.0, nq + 1))


def quad_on(fun, edges):
    """Composite 8-point Gauss with the integration cells given by `edges`.

    Needed whenever the integrand is only piecewise smooth: the derivative of a
    P1 function jumps at the nodes, and a cell that straddles a node would make
    the integral inexact no matter how many Gauss points are used.
    """
    edges = np.asarray(edges, float)
    a = edges[:-1][:, None]
    b = edges[1:][:, None]
    xm, xr = 0.5 * (a + b), 0.5 * (b - a)
    x = (xm + xr * _GP[None, :]).ravel()
    w = (xr * _GW[None, :]).ravel()
    return float(np.dot(w, fun(x)))


def exact_u(x):
    return np.sin(np.pi * x)


def exact_du(x):
    return np.pi * np.cos(np.pi * x)


def load(x):
    return np.pi ** 2 * np.sin(np.pi * x)


def uh(x, h, c):
    """Piecewise-linear hat interpolant and its (piecewise constant) derivative."""
    n = len(c)
    xi = x / h
    i = np.clip(np.floor(xi).astype(int), 0, n)
    s = xi - i
    left = np.where(i >= 1, c[np.clip(i - 1, 0, n - 1)], 0.0)
    right = np.where(i <= n - 1, c[np.clip(i, 0, n - 1)], 0.0)
    return (1.0 - s) * left + s * right, (right - left) / h


def assemble(n):
    """Uniform mesh, n interior nodes, h = 1/(n+1); returns h, K, f."""
    h = 1.0 / (n + 1)
    nodes = np.linspace(0.0, 1.0, n + 2)
    Ke = np.array([[1.0, -1.0], [-1.0, 1.0]]) / h
    K = np.zeros((n, n))
    for e in range(n + 1):
        idx = [e - 1, e]
        for p in range(2):
            for q in range(2):
                if 0 <= idx[p] < n and 0 <= idx[q] < n:
                    K[idx[p], idx[q]] += Ke[p, q]
    f = np.zeros(n)
    for e in range(n + 1):
        a, b = nodes[e], nodes[e + 1]
        xm, xr = 0.5 * (a + b), 0.5 * (b - a)
        x = xm + xr * _GP
        w = xr * _GW
        s = (x - a) / h
        fe = np.array([np.dot(w, load(x) * (1.0 - s)),
                       np.dot(w, load(x) * s)])
        for p in range(2):
            if 0 <= e - 1 + p < n:
                f[e - 1 + p] += fe[p]
    return h, K, f


def fe_solution(n):
    h, K, f = assemble(n)
    return h, K, f, np.linalg.solve(K, f)


# ================================================================ 1. Ritz
def part1():
    print("\n# 1  Ritz: J'(c) = K c - f = 0, and J is minimal there")
    h, K, f, c = fe_solution(7)
    check_true("K is symmetric", np.allclose(K, K.T, atol=0),
               "max|K-K^T| = %.1e" % np.abs(K - K.T).max())
    check_true("K is positive definite", np.all(np.linalg.eigvalsh(K) > 0),
               "min eig = %.4f" % np.linalg.eigvalsh(K).min())
    check("residual |K c - f|", np.linalg.norm(K @ c - f), 0.0, 1e-12)

    def J(cc):
        return 0.5 * cc @ K @ cc - f @ cc

    rng = np.random.default_rng(20261113)
    worst = min(J(c + 1e-3 * rng.standard_normal(len(c))) - J(c) for _ in range(200))
    check_true("J(c+dc) > J(c) for 200 random dc", worst > 0,
               "min increase = %.3e" % worst)
    check("J at the Ritz point = -A(uN,uN)/2", J(c), -0.5 * c @ K @ c, 1e-12)


# ============================================================ 2. Galerkin
def part2():
    print("\n# 2  Galerkin: A(uN,vN) = F(vN) for every vN in S_N")
    for n in (3, 7, 31):
        h, K, f, c = fe_solution(n)
        rng = np.random.default_rng(99 + n)
        nodes = np.linspace(0.0, 1.0, n + 2)
        for trial in range(3):
            w = rng.standard_normal(n)
            lhs = quad_on(lambda x: uh(x, h, c)[1] * uh(x, h, w)[1], nodes)
            rhs = quad_on(lambda x: load(x) * uh(x, h, w)[0], nodes)
            check("n=%-2d trial %d  |A(uN,vN)-F(vN)| / (|A|+|F|)"
                  % (n, trial), abs(lhs - rhs) / (abs(lhs) + abs(rhs)), 0.0, 1e-10)


# ============================================== 3. energy identity / lower bound
def part3():
    print("\n# 3  A(u,u) - A(uN,uN) = ||u-uN||_A^2  (identity, hence a lower bound)")
    for n in (3, 7, 31):
        h, K, f, c = fe_solution(n)
        A_h = c @ K @ c
        nodes = np.linspace(0.0, 1.0, n + 2)
        err2 = quad_on(lambda x: (exact_du(x) - uh(x, h, c)[1]) ** 2, nodes)
        check("identity  n=%-2d" % n, A_EXACT - A_h, err2, 1e-9)
        check_true("A(uN,uN) < A(u,u)  n=%-2d" % n, A_h < A_EXACT,
                   "A(uN,uN) = %.8f  <  %.8f" % (A_h, A_EXACT))


# ============================================ 4. orthogonality / best approximation
def part4():
    print("\n# 4  orthogonality A(u-uN, phi_i) = 0  and best approximation")
    n = 15
    h, K, f, c = fe_solution(n)
    nodes = np.linspace(0.0, 1.0, n + 2)

    worst = 0.0
    for i in range(1, n + 1):
        e_i = np.zeros(n)
        e_i[i - 1] = 1.0
        res = quad_on(lambda x: (exact_du(x) - uh(x, h, c)[1]) * uh(x, h, e_i)[1],
                      nodes)
        worst = max(worst, abs(res))
    check("max |A(u-uN, phi_i)|  n=%d" % n, worst, 0.0, 1e-9)

    def energy_err(cc):
        return quad_on(lambda x: (exact_du(x) - uh(x, h, cc)[1]) ** 2, nodes)

    e0 = energy_err(c)
    rng = np.random.default_rng(7)
    worst = min(energy_err(c + 1e-3 * rng.standard_normal(n)) - e0 for _ in range(100))
    check_true("100 random competitors are farther in the energy norm", worst > 0,
               "min increase = %.3e" % worst)


# ================================================================ 5. rates
def part5():
    print("\n# 5  convergence table quoted in section 3.3")
    print("    n      h         A(uN,uN)       eA          ratio    eL2         ratio")
    rows, prev = [], None
    for n in (3, 7, 15, 31, 63, 127):
        h, K, f, c = fe_solution(n)
        A_h = c @ K @ c
        eA = np.sqrt(A_EXACT - A_h)
        nodes = np.linspace(0.0, 1.0, n + 2)
        eL2 = np.sqrt(quad_on(lambda x: (exact_u(x) - uh(x, h, c)[0]) ** 2, nodes))
        rA = "" if prev is None else "%.3f" % (prev[0] / eA)
        rL = "" if prev is None else "%.3f" % (prev[1] / eL2)
        print("    %-6d %-9.6f %-14.9f %-11.3e %-8s %-11.3e %-8s"
              % (n, h, A_h, eA, rA, eL2, rL))
        rows.append((n, h, A_h, eA, eL2))
        prev = (eA, eL2)

    for a, b in zip(rows[:-1], rows[1:]):
        check("energy-norm rate n=%d->%d (want ~2)" % (a[0], b[0]),
              a[3] / b[3], 2.0, 0.12)
        check("L2 rate n=%d->%d (want ~4)" % (a[0], b[0]),
              a[4] / b[4], 4.0, 0.45)
    return rows


# ==================================================== 6. nonsymmetric 2x2
def part6():
    print("\n# 6  nonsymmetric A: Ritz (symmetrised) and Galerkin differ")
    A = np.array([[2.0, 1.0], [0.0, 3.0]])
    F = np.array([1.0, 0.0])
    c_gal = np.linalg.solve(A, F)
    c_ritz = np.linalg.solve(0.5 * (A + A.T), F)
    print("    Galerkin c = %s" % np.array2string(c_gal, precision=6))
    print("    Ritz     c = %s" % np.array2string(c_ritz, precision=6))
    check("Galerkin c1", c_gal[0], 0.5, 1e-12)
    check("Galerkin c2", c_gal[1], 0.0, 1e-12)
    check("Ritz c1", c_ritz[0], 12.0 / 23.0, 1e-12)
    check("Ritz c2", c_ritz[1], -2.0 / 23.0, 1e-12)
    check("||c_gal - c_ritz||_2", np.linalg.norm(c_gal - c_ritz), 0.08963273, 1e-7)


# ==================================================== 7. sparsity / cost
def part7():
    print("\n# 7  sparsity of the hat-function stiffness matrix")
    for n in (31, 99, 999):
        nnz = 3 * n - 2
        print("    n=%-4d nnz=%-5d density=%.4f%%"
              % (n, nnz, 100.0 * nnz / n ** 2))
    check("nnz  n=999", 3 * 999 - 2, 2995, 0)
    check("density  n=99  (%)", 100.0 * (3 * 99 - 2) / 99 ** 2, 3.0, 0.05)
    check("density  n=999 (%)", 100.0 * (3 * 999 - 2) / 999 ** 2, 0.3, 0.005)


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
    part7()
    print("\nall variational-note numbers verified")
