#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""校验 FRIC 笔记（M6-09）中的全部数值结论。

运行：
    py -3 tools/verify_fric_notes.py

末尾打印 ALL CHECKS PASSED；任一检查失败立即退出。
"""

import numpy as np
import fric_lib as fl


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
    print("\n[1] 算例一：三个增量、三个方向")
    rows = fl.example1_steps()
    check("步 3 滑移 x", rows[0]["dgam"][0], 0.1, 1e-12)
    check("步 3 滑移 y", rows[0]["dgam"][1], 0.0, 1e-12)
    approx("步 4 滑移 x", rows[1]["dgam"][0], 0.1 / np.sqrt(2.0), 1e-12)
    approx("步 4 滑移 y", rows[1]["dgam"][1], 0.1 / np.sqrt(2.0), 1e-12)
    check("步 5 滑移 x", rows[2]["dgam"][0], 0.0, 1e-12)
    check("步 5 滑移 y", rows[2]["dgam"][1], 0.1, 1e-12)
    approx("终点 x", rows[2]["cum"][0], 0.17071, 1e-4)
    approx("终点 y", rows[2]["cum"][1], 0.17071, 1e-4)
    check("路径长", rows[2]["path"], 0.3, 1e-12)
    s = fl.example1_summary()
    approx("弦长", s["chord"], 0.2414213562373095, 1e-12)
    check("方向角 (deg)", s["angle"], 45.0, 1e-10)
    approx("弦/路径", s["chord"] / s["path"], 0.804738, 1e-6)
    for r in rows:
        approx("每步增量模长", np.linalg.norm(r["dgam"]), 0.1, 1e-12)
    print("  终点 = (%.8f, %.8f), 路径 = %.6f, 弦 = %.6f"
          % (s["final"][0], s["final"][1], s["path"], s["chord"]))


def part2():
    print("\n[2] 算例二：温度依赖 k 与增量采样")
    check("k(0)", fl.example2_k(0.0), 0.001, 1e-15)
    check("k(1)", fl.example2_k(1.0), 0.002, 1e-15)
    check("连续解", fl.example2_continuous(), 0.15, 1e-12)
    expected_end = {1: 0.2, 2: 0.175, 5: 0.16, 10: 0.155,
                    20: 0.1525, 100: 0.1505}
    for n, target in expected_end.items():
        check("末采样 n=%d" % n, fl.example2_end_slip(n), target, 1e-12)
        check("末采样闭式 n=%d" % n,
              fl.example2_end_slip_closed(n), target, 1e-12)
        check("中点采样 n=%d" % n, fl.example2_mid_slip(n), 0.15, 1e-12)
        check("中点采样闭式 n=%d" % n,
              fl.example2_mid_slip_closed(n), 0.15, 1e-12)
    check("末采样 n=10 与连续解之差",
          fl.example2_end_slip(10) - fl.example2_continuous(), 0.005, 1e-12)
    approx("末采样 n=10 相对误差",
          (fl.example2_end_slip(10) - 0.15) / 0.15, 1.0 / 30.0, 1e-12)
    ratio = ((fl.example2_end_slip(10) - 0.15)
             / (fl.example2_end_slip(20) - 0.15))
    check("误差比 e(10)/e(20)", ratio, 2.0, 1e-10)
    for n in (1, 2, 5, 10, 20, 100):
        check("一阶收敛 e(n)*n", (fl.example2_end_slip(n) - 0.15) * n,
              0.05, 1e-12)


def part3():
    print("\n[3] 切线：ddtddg 与 ddtddt")
    tau, ddtddg, ddtddt = fl.fric_tangents(3e-4, 0.2, theta_bar=150.0)
    check("tau", tau, 0.6, 1e-12)
    check("ddtddg", ddtddg, 2000.0, 1e-10)
    check("ddtddt", ddtddt, -0.0012, 1e-15)

    h = 1e-8
    tau_p = fl.fric_tangents(3e-4 + h, 0.2, theta_bar=150.0)[0]
    tau_m = fl.fric_tangents(3e-4 - h, 0.2, theta_bar=150.0)[0]
    approx("ddtddg 有限差分", (tau_p - tau_m) / (2 * h), 2000.0, 1e-8)

    def tau_of(theta1, theta2):
        k = 0.001 + 1e-5 * 0.5 * (theta1 + theta2)
        return 3e-4 / (k * 0.2)

    h2 = 1e-3
    fd1 = (tau_of(150.0 + h2, 150.0) - tau_of(150.0 - h2, 150.0)) / (2 * h2)
    fd2 = (tau_of(150.0, 150.0 + h2) - tau_of(150.0, 150.0 - h2)) / (2 * h2)
    approx("ddtdtheta1 有限差分", fd1, -0.0012, 1e-8)
    approx("ddtdtheta2 有限差分", fd2, -0.0012, 1e-8)
    check("两侧温度导数相等", fd1 - fd2, 0.0, 1e-12)


def part4():
    print("\n[4] 摩擦热账本")
    check("t=0 热流", fl.friction_heat(0.001, 100.0), 5.0, 1e-12)
    check("t=1 热流", fl.friction_heat(0.002, 100.0), 10.0, 1e-12)
    check("连续累计热", fl.example2_heat_continuous(), 7.5, 1e-12)
    check("末采样 n=10 累计热", fl.example2_heat_end(10), 7.75, 1e-12)
    check("中点采样 n=10 累计热", fl.example2_heat_mid(10), 7.5, 1e-12)
    check("中点采样 n=100 累计热", fl.example2_heat_mid(100), 7.5, 1e-12)
    work_c = 100.0 * fl.example2_continuous()
    work_e = 100.0 * fl.example2_end_slip(10)
    check("连续摩擦功", work_c, 15.0, 1e-12)
    check("末采样 n=10 摩擦功", work_e, 15.5, 1e-12)
    approx("热/功 = eta", fl.example2_heat_continuous() / work_c, 0.5, 1e-12)
    check("连续每个表面", fl.example2_heat_continuous() / 2.0, 3.75, 1e-12)
    check("末采样每个表面", fl.example2_heat_end(10) / 2.0, 3.875, 1e-12)
    approx("每个表面/功", (fl.example2_heat_end(10) / 2.0) / work_e,
           0.25, 1e-12)
    print("  连续：功 15.0，热 7.5，每面 3.75；"
          "末采样 n=10：功 15.5，热 7.75，每面 3.875")


def part5():
    print("\n[5] 库仑对照与量纲")
    check("库仑上限 mu*p", fl.coulomb_limit(0.6, 100.0), 60.0, 1e-12)
    assert fl.coulomb_regime(40.0, 0.6, 100.0) == "stick"
    assert fl.coulomb_regime(60.0, 0.6, 100.0) == "slip"
    assert fl.coulomb_regime(80.0, 0.6, 100.0) == "slip"
    print("PASS  库仑判据 tau<60 粘着，tau>=60 滑移")
    check("线性黏性 tau=40 的滑移速度", fl.slip_rate(40.0, 0.001), 0.04,
          1e-15)
    check("线性黏性 tau=0 的滑移速度", fl.slip_rate(0.0, 0.001), 0.0, 1e-15)
    check("线性黏性 tau=1e-6 的滑移速度",
          fl.slip_rate(1e-6, 0.001), 1e-9, 1e-20)
    print("PASS  线性黏性没有粘着门槛：任意非零 tau 都产生非零滑移速度")

    # 单位：同一 k、同一 tau，长度单位换成 mm 只改变 dgam 的单位
    dgam_m = fl.slip_increment(100.0, 0.001, 1.0)
    check("SI 长度单位下的 dgam (m)", dgam_m, 0.1, 1e-12)
    check("换算成 mm", dgam_m * 1000.0, 100.0, 1e-12)
    # 若把 k 当作 1/(Pa*s)，得到的是无量纲剪应变增量；
    # 乘以参考长度 1 mm 才回到 Abaqus 的长度口径
    strain_inc = 0.001 * 100.0 * 1.0
    check("无量纲剪应变增量", strain_inc, 0.1, 1e-12)
    check("乘参考长度 1 mm 后的滑移 (mm)", strain_inc * 1.0, 0.1, 1e-12)

    # 力/应力口径：面积缩小 100 倍，tau 与 dgam 都放大 100 倍
    tau_unit = 100.0 / 1.0
    tau_small = 100.0 / 0.01
    check("单位面积 tau", tau_unit, 100.0, 1e-12)
    check("0.01 m^2 上的 tau", tau_small, 10000.0, 1e-12)
    check("面积缩小后的 dgam",
          fl.slip_increment(tau_small, 0.001, 1.0), 10.0, 1e-12)
    print("PASS  力与应力的数值只在单位面积时相同；面积 0.01 m^2 时差 100 倍")


def main():
    part1()
    part2()
    part3()
    part4()
    part5()
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
