#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recompute every number quoted in docs/06-numerics/element-types.md.

Usage (from the repository root):
    py -3 tools/verify_element_notes.py

单位：mm、MPa、N。脚本逐项打印，任一断言失败即以非零码退出。
"""

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import element_lib as el  # noqa: E402


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


def q4_stiff(coords):
    return el.isoparametric_stiffness(coords, el.q4_shape, el.gauss_2d(2))


def q4_reduced(coords):
    return el.isoparametric_stiffness(coords, el.q4_shape, el.gauss_2d(1))


def q8_stiff(coords):
    return el.isoparametric_stiffness(coords, el.q8_shape, el.gauss_2d(3))


# ============================================================ 1. 形函数与基准
def part1():
    print("\n# 1  形函数、弹性矩阵与梁理论基准")
    for xi, eta in [(0.37, -0.21), (-0.6, 0.8), (0.0, 0.0)]:
        for name, shape in [("Q4", el.q4_shape), ("Q8", el.q8_shape)]:
            N, dN = shape(xi, eta)
            check("%s sum(N) at (%.2f,%.2f)" % (name, xi, eta),
                  N.sum(), 1.0, 1e-14)
            check("%s sum(dN/dxi)" % name, dN[:, 0].sum(), 0.0, 1e-13)
            check("%s sum(dN/deta)" % name, dN[:, 1].sum(), 0.0, 1e-13)

    D = el.plane_stress_D()
    check("D symmetric", D[0, 1] - D[1, 0], 0.0, 1e-9)
    check("D shear term E/(2(1+nu))", D[2, 2],
          el.E_STEEL / (2 * (1 + el.NU_STEEL)), 1e-8)
    approx("Euler-Bernoulli tip deflection (mm)", el.beam_tip_exact(),
           1.9047619047619047, 1e-12)


# ============================================================ 2. 悬臂梁精度
def part2():
    print("\n# 2  悬臂梁：Q4 / Q8 / T3 的端部挠度")
    exact = el.beam_tip_exact()

    nodes, elems = el.mesh_q4(10, 1)
    check("Q4 10x1 node count", len(nodes), 22, 0)
    check("Q4 10x1 element count", len(elems), 10, 0)
    approx("Q4 10x1 / exact", el.tip_deflection(nodes, elems, q4_stiff) / exact,
           0.67889, 0.001)

    for ny, target in [(2, 0.71109), (4, 0.72015), (8, 0.72247)]:
        nodes, elems = el.mesh_q4(10, ny)
        approx("Q4 10x%d / exact" % ny,
               el.tip_deflection(nodes, elems, q4_stiff) / exact, target, 0.001)

    nodes, elems = el.mesh_q4(80, 8)
    check("Q4 80x8 node count", len(nodes), 729, 0)
    check("Q4 80x8 element count", len(elems), 640, 0)
    approx("Q4 80x8 / exact (640 elements)",
           el.tip_deflection(nodes, elems, q4_stiff) / exact, 0.99796, 0.0005)

    nodes, elems = el.mesh_q8(10, 1)
    check("Q8 10x1 node count", len(nodes), 53, 0)
    check("Q8 10x1 element count", len(elems), 10, 0)
    approx("Q8 10x1 / exact (one layer)",
           el.tip_deflection(nodes, elems, q8_stiff) / exact, 0.99826, 0.0005)

    nodes, elems = el.mesh_t3(10, 1)
    check("T3 10x1 element count", len(elems), 20, 0)
    approx("T3 10x1 / exact",
           el.tip_deflection(nodes, elems, el.cst_stiffness) / exact,
           0.23118, 0.001)


# ============================================================ 3. 单层天花板
def part3():
    print("\n# 3  只沿长度加密：单层 Q4 的天花板")
    exact = el.beam_tip_exact()
    for nx, target in [(20, 0.84276), (40, 0.89688), (160, 0.91525)]:
        nodes, elems = el.mesh_q4(nx, 1)
        approx("Q4 %dx1 / exact" % nx,
               el.tip_deflection(nodes, elems, q4_stiff) / exact, target, 0.001)


# ============================================================ 4. 剪切锁死
def part4():
    print("\n# 4  纯弯曲下的寄生剪切（剪切锁死）")
    gamma = el.parasitic_shear(kappa=1.0, a=10.0)
    approx("analytic gamma at a/h=10", gamma, 2.886751345948129, 1e-10)

    strains = el.pure_bending_strains(a=10.0, h=1.0, kappa=1.0)
    approx("eps_x at lower Gauss point", strains[0][0], 0.2886751346, 1e-9)
    approx("eps_x at upper Gauss point", strains[1][0], -0.2886751346, 1e-9)
    check("gamma sign flips between left and right columns",
          strains[0][1] + strains[2][1], 0.0, 1e-12)
    approx("|gamma| equals kappa a/(2 sqrt3)", abs(strains[0][1]),
           2.886751345948129, 1e-10)
    approx("|gamma / eps_x| at lower point", abs(strains[0][1] / strains[0][0]),
           10.0, 1e-8)

    check("Q8 pure-bending shear is zero",
          el.q8_pure_bending_shear(a=10.0, h=1.0), 0.0, 1e-10)


# ============================================================ 5. 特征值
def part5():
    print("\n# 5  单刚特征值：刚体模态与沙漏模态")
    cases = [
        ("Q4 2x2", el.q4_shape, el.gauss_2d(2), 3),
        ("Q4 1x1", el.q4_shape, el.gauss_2d(1), 5),
        ("Q8 3x3", el.q8_shape, el.gauss_2d(3), 3),
        ("Q8 2x2", el.q8_shape, el.gauss_2d(2), 4),
    ]
    stored = {}
    for name, shape, gauss, nzero in cases:
        eigs = el.element_eigenvalues(shape, gauss)
        stored[name] = eigs
        check("%s zero-mode count" % name,
              el.zero_mode_count(eigs), nzero, 0)

    approx("Q4 2x2 smallest positive eigenvalue",
           stored["Q4 2x2"][3], 1.0385e5, 2e-4)
    approx("Q4 1x1 smallest positive eigenvalue",
           stored["Q4 1x1"][5], 1.6154e5, 2e-4)
    approx("Q8 3x3 smallest positive eigenvalue",
           stored["Q8 3x3"][3], 3.5291e4, 2e-4)
    approx("Q8 2x2 smallest positive eigenvalue",
           stored["Q8 2x2"][4], 6.3346e4, 2e-4)

    modes = el.hourglass_modes()
    check("hourglass mode count", len(modes), 2, 0)
    K = el.isoparametric_stiffness(
        np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
        el.q4_shape, el.gauss_2d(1), thickness=1.0)
    check("K_reduced @ hourglass mode = 0",
          np.max(np.abs(K @ modes.T)), 0.0, 1e-8)
    check("hourglass modes normalized to max 1",
          np.max(np.abs(modes)), 1.0, 1e-12)

    # 沙漏模态必须与刚体模态正交：平移 (1,0)、(0,1)，转动 (-y, x)
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    R = np.zeros((8, 3))
    for k, (x, y) in enumerate(coords):
        R[2 * k, 0] = 1.0
        R[2 * k + 1, 1] = 1.0
        R[2 * k, 2] = -y
        R[2 * k + 1, 2] = x
    check("hourglass orthogonal to rigid modes",
          np.max(np.abs(R.T @ modes.T)), 0.0, 1e-10)


# ============================================================ 6. 工程账
def part6():
    print("\n# 6  挡泥板：实体网格 vs 壳网格")
    acc = el.dof_accounting()
    s, sh = acc["solid"], acc["shell"]
    check("solid element count", s["elements"], 240000, 0)
    check("solid node count", s["nodes"], 303505, 0)
    check("solid DOF (3 per node)", s["dof"], 910515, 0)
    check("shell element count", sh["elements"], 2400, 0)
    check("shell node count", sh["nodes"], 2541, 0)
    check("shell DOF (6 per node)", sh["dof"], 15246, 0)
    approx("solid / shell element ratio",
           s["elements"] / sh["elements"], 100.0, 1e-12)
    approx("solid / shell DOF ratio", s["dof"] / sh["dof"], 59.72104, 1e-4)

    pcr = el.spoke_buckling()
    approx("spoke Euler buckling load (N)", pcr, 19.35588, 1e-4)
    approx("spoke Euler buckling load (kg)", pcr / 9.80665, 1.97375, 1e-4)
    approx("spoke tensile stress at 100 kg (MPa)",
           el.spoke_tension_stress(100.0, 2.0), 312.155, 1e-4)


# ============================================================ 7. 二维降维
def part7():
    print("\n# 7  平面应力与平面应变的 sigma_z")
    approx("plane strain sigma_z (sx=sy=-10 MPa)",
           el.plane_strain_sigma_z(-10.0, -10.0), -6.0, 1e-12)
    approx("plane stress eps_z (sx=sy=-10 MPa)",
           el.plane_stress_eps_z(-10.0, -10.0), 2.857142857e-5, 1e-10)
    check("plane stress sigma_z is zero by definition",
          el.plane_strain_sigma_z(0.0, 0.0), 0.0, 1e-15)


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
    part7()
    print("\nALL CHECKS PASSED")
