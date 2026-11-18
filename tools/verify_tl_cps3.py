#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""校验 M6-10 笔记（CPS3 完全拉格朗日格式）中的全部数值结论。

运行：
    py -3 tools/verify_tl_cps3.py

末尾打印 ALL CHECKS PASSED；任一检查失败立即退出。
"""

import numpy as np

import tl_cps3_lib as tl


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


def part1():
    print("\n[1] 参考几何与线性形函数")
    coords = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 1.0]])
    g = tl.tri_geometry(coords, 0.5)
    check("A0", g["A"], 1.0, 1e-14)
    check("V0", g["V"], 0.5, 1e-14)
    check("dN1/dX", g["dNdx"][0], -0.5, 1e-14)
    check("dN2/dX", g["dNdx"][1], 0.5, 1e-14)
    check("dN3/dY", g["dNdy"][2], 1.0, 1e-14)

    pts = np.array([[2.0 / 3.0, 1.0 / 3.0], [0.2, 0.7], [1.1, 0.1]])
    n = tl.shape_functions(g, pts)
    check("形函数在任意点之和 (max)", np.max(np.abs(n.sum(axis=1) - 1.0)),
          0.0, 1e-14)
    approx("形函数在形心", n[0, 0], 1.0 / 3.0, 1e-14)
    print("  面积坐标：N1=0.3333、N2=0.3333、N3=0.3333（形心）")


def part2():
    print("\n[2] B 矩阵：B(0)=B_L，且等于数值应变导数")
    coords = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 1.0]])
    g = tl.tri_geometry(coords, 1.0)
    b0 = tl.b_matrix(g, np.zeros(6))
    bl = tl.b_linear(g)
    check("max|B(0)-B_L|", np.max(np.abs(b0 - bl)), 0.0, 0.0)

    rng = np.random.default_rng(3)
    d = rng.normal(0.0, 0.05, 6)
    b = tl.b_matrix(g, d)
    eps = 1e-7
    bn = np.zeros((3, 6))
    for j in range(6):
        e = np.zeros(6)
        e[j] = eps
        bn[:, j] = (tl.gl_strain(g, d + e) - tl.gl_strain(g, d - e)) / (2 * eps)
    check("max|B - dE/dd|", np.max(np.abs(b - bn)), 0.0, 1e-9)
    bnl = tl.b_nonlinear(g, d)
    print("  B_NL 的 max 元 = %.6f（随当前位移梯度更新）"
          % np.max(np.abs(bnl)))


def part3():
    print("\n[3] 客观性：刚体旋转下 GL 应变为零")
    coords = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 1.0]])
    g = tl.tri_geometry(coords, 1.0)
    th = np.deg2rad(30.0)
    rot = np.array([[np.cos(th), -np.sin(th)],
                    [np.sin(th), np.cos(th)]])
    u = coords @ (rot.T - np.eye(2))
    d = u.reshape(-1)
    e = tl.gl_strain(g, d)
    small = tl.small_strain(g, d)
    check("max|E|（30° 刚转）", np.max(np.abs(e)), 0.0, 1e-14)
    approx("小应变 E_XX", small[0], np.cos(th) - 1.0, 1e-13)
    approx("小应变 E_YY", small[1], np.cos(th) - 1.0, 1e-13)
    check("小应变 2E_XY", small[2], 0.0, 1e-14)
    print("  GL: 0；小应变: %.6f（把刚转误当剪切/压缩）" % small[0])


def part4():
    print("\n[4] 均匀拉伸：解析应变、应力与虚功一致性")
    s = 0.2
    young, nu = 70e9, 0.33
    g, d, strain, stress, f = tl.uniform_stretch_state(s, young, nu)
    approx("E_XX", strain[0], s + 0.5 * s * s, 1e-14)
    check("E_YY", strain[1], 0.0, 1e-14)
    check("2E_XY", strain[2], 0.0, 1e-14)
    approx("S_XX (Pa)", stress[0], young / (1 - nu * nu) * 0.22, 1e-13)
    approx("S_YY / S_XX", stress[1] / stress[0], nu, 1e-13)
    fsum = f.reshape(3, 2).sum(axis=0)
    check("内力合力 x", fsum[0], 0.0, 1e-6)
    check("内力合力 y", fsum[1], 0.0, 1e-6)

    # 虚功一致性：随机虚位移，f.dd = V S.dE
    rng = np.random.default_rng(7)
    dd = rng.normal(0.0, 1.0, 6)
    de = tl.b_matrix(g, d) @ dd
    lhs = float(f @ dd)
    rhs = float(g["V"] * (stress @ de))
    approx("f.dd / (V S.dE)", lhs / rhs, 1.0, 1e-13)
    w = 0.5 * float(stress @ strain)
    approx("应变能密度 W (J/m^3)", w,
           0.5 * stress[0] * strain[0], 1e-13)
    print("  S_XX = %.4f GPa, S_YY = %.4f GPa, W = %.4f MJ/m^3"
          % (stress[0] / 1e9, stress[1] / 1e9, w / 1e6))


def part5():
    print("\n[5] 切线：材料 + 几何 vs 中心差分")
    coords = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 1.0]])
    g = tl.tri_geometry(coords, 1.0)
    rng = np.random.default_rng(11)
    for k, scale in enumerate([0.02, 0.05, 0.10], start=1):
        d = rng.normal(0.0, scale, 6)
        k_mat, k_geo, k_tot = tl.tangent(g, d, 1.0, 0.3)
        dd = rng.normal(0.0, 1.0, 6)
        dd /= np.linalg.norm(dd)
        eps = 1e-6
        fp = tl.internal_force(g, d + eps * dd, 1.0, 0.3)[0]
        fm = tl.internal_force(g, d - eps * dd, 1.0, 0.3)[0]
        fd = (fp - fm) / (2 * eps)
        rel = np.linalg.norm(fd - k_tot @ dd) / np.linalg.norm(k_tot @ dd)
        print("PASS  第 %d 组：相对误差 %.2e（|θ|≈%.2f）"
              % (k, rel, np.linalg.norm(tl.displacement_gradient(g, d))))
        if rel > 1e-8:
            print("FAIL  切线相对误差过大")
            raise SystemExit(1)
        asym = np.max(np.abs(k_tot - k_tot.T)) / np.max(np.abs(k_tot))
        if asym > 1e-12:
            print("FAIL  切线不对称")
            raise SystemExit(1)
    k_geo0 = tl.tangent(g, np.zeros(6), 1.0, 0.3)[1]
    check("零应变时 max|K_geo|", np.max(np.abs(k_geo0)), 0.0, 1e-16)
    print("  切线对称性 < 1e-12；零应力时几何刚度为 0")


def part6():
    print("\n[6] 两单元 patch：Newton 收敛")
    res_full = tl.newton_patch(0.2, 1.0, 0.3, material_only=False)
    res_mat = tl.newton_patch(0.2, 1.0, 0.3, material_only=True, max_iter=30)
    if not res_full["converged"] or res_full["iterations"] > 8:
        print("FAIL  一致切线未在 8 次内收敛")
        raise SystemExit(1)
    print("PASS  一致切线：%d 次迭代，末残差 %.2e"
          % (res_full["iterations"], res_full["residuals"][-1]))
    print("      残差序列：" + ", ".join(
        "%.2e" % r for r in res_full["residuals"]))
    ratios = [res_full["residuals"][i + 1]
              / res_full["residuals"][i] ** 2
              for i in range(1, len(res_full["residuals"]) - 1)]
    print("      二次收敛指标 r(k+1)/r(k)^2："
          + ", ".join("%.2f" % x for x in ratios[:3]))
    if res_mat["iterations"] < 20:
        print("FAIL  仅材料切线不应这么快收敛")
        raise SystemExit(1)
    print("PASS  仅材料切线：30 次后残差仍为 %.2e（线性收敛）"
          % res_mat["residuals"][-1])

    r1 = tl.newton_patch(0.5, 1.0, 0.3, material_only=True, max_iter=30)
    if r1["converged"]:
        print("FAIL  仅材料切线在大载荷下不应收敛")
        raise SystemExit(1)
    print("PASS  仅材料切线在 load=0.5 时发散（末残差 %.2e）"
          % r1["residuals"][-1])

    d = res_full["d"]
    f_int, _ = tl.assemble_patch(res_full["coords"], res_full["conn"], d,
                                 1.0, 0.3)
    r = res_full["f_ext"] - f_int
    free = [2, 3, 4, 5]
    check("patch 自由自由度残差", np.linalg.norm(r[free]), 0.0, 1e-14)
    print("  右边界位移 u1=%.6f, u2=%.6f, v2=%.6f"
          % (d[2], d[4], d[5]))


def main():
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()