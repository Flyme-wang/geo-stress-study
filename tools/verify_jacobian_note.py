#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recompute every number quoted in docs/06-numerics/jacobian-and-mapping.md.

Usage (from the repository root):
    py -3 tools/verify_jacobian_note.py

The note quotes det J at four working points, four area-ratio tables, a
first-order error table, the O(h^2) behaviour near a singular point, and the
det J = l1*l2*sin(theta2) of a planar two-link arm. Every one of them is
recomputed here; the script exits non-zero on the first mismatch.
"""

import sys

import numpy as np


def f(x):
    """Running example of the note: y1 = x1^2 + x2, y2 = x1*x2."""
    return np.array([x[0] ** 2 + x[1], x[0] * x[1]], float)


def jac(x):
    return np.array([[2.0 * x[0], 1.0], [x[1], x[0]]])


def check(label, value, target, tol):
    ok = abs(value - target) <= tol
    print(("PASS  " if ok else "FAIL  ") + label
          + "  value=%.6f  target=%.6f  tol=%.1e" % (value, target, tol))
    if not ok:
        raise SystemExit(1)


def area_ratio(x0, h):
    """Area of the image of an h x h square, divided by h^2 (shoelace)."""
    sq = np.array([[0, 0], [h, 0], [h, h], [0, h]], float)
    pts = np.array([f(x0 + s) for s in sq])
    area = 0.5 * abs(np.dot(pts[:, 0], np.roll(pts[:, 1], -1))
                     - np.dot(pts[:, 1], np.roll(pts[:, 0], -1)))
    return area / h ** 2


# ============================================================ 1. det J
def part1():
    print("\n# 1  det J at the four working points")
    expected = [((2.0, 1.0), 7.0), ((0.5, 0.0), 0.5),
                ((0.0, 1.0), -1.0), ((1.0, 2.0), 0.0)]
    for p, target in expected:
        det = np.linalg.det(jac(np.array(p, float)))
        check("det J at %s" % (p,), det, target, 1e-12)


# =================================================== 2. area ratios
def part2():
    print("\n# 2  area ratio -> |det J|, table in section 3.2")
    table = {
        (2.0, 1.0): [7.355, 7.0351, 7.0035, 7.00035],
        (0.5, 0.0): [0.555, 0.50505, 0.5005, 0.50005],
        (0.0, 1.0): [1.045, 1.00495, 1.0005, 1.00005],
        (1.0, 2.0): [0.155, 0.01505, 0.0015005, 0.000149925],
    }
    hs = [1e-1, 1e-2, 1e-3, 1e-4]
    for p, row in table.items():
        x0 = np.array(p, float)
        for h, target in zip(hs, row):
            check("area ratio %s h=%.0e" % (p, h), area_ratio(x0, h),
                  target, 5e-7 * max(1.0, abs(target)) / max(1e-4, h)
                  if target else 1e-9)

    # the non-singular columns converge to |det| with error ~ C*h
    for p in [(2.0, 1.0), (0.5, 0.0), (0.0, 1.0)]:
        x0 = np.array(p, float)
        det = abs(np.linalg.det(jac(x0)))
        r1, r2 = area_ratio(x0, 1e-2), area_ratio(x0, 1e-3)
        check("area ratio error falls ~10x  %s" % (p,),
              (r1 - det) / (r2 - det), 10.0, 1.0)

    # the singular point: area is O(h^3), i.e. ratio/h is bounded and ~1.5
    x0 = np.array([1.0, 2.0])
    ratio = area_ratio(x0, 1e-3)
    check("singular point ratio/h", ratio / 1e-3, 1.5, 0.02)


# ============================================== 3. first-order error
def part3():
    print("\n# 3  first-order error along (h, 0) at (2, 1)")
    x0 = np.array([2.0, 1.0])
    table = [(1e-1, 2.42536e-2), (1e-2, 2.42536e-3),
             (1e-3, 2.42536e-4), (1e-4, 2.42536e-5)]
    for h, target in table:
        d = np.array([h, 0.0])
        # exact increment written analytically: y1 = (2+h)^2+1-5, y2 = (2+h)*1-2.
        # Evaluating f() directly loses digits to cancellation at h = 1e-4.
        exact = np.array([4.0 * h + h * h, h])
        approx = jac(x0) @ d
        rel = np.linalg.norm(exact - approx) / np.linalg.norm(approx)
        check("relative error h=%.0e (= h/sqrt(17))" % h,
              rel, h / np.sqrt(17), 1e-12)
        check("table entry h=%.0e" % h, rel, target, 1e-6)


# ======================================== 4. singular point, O(h^2)
def part4():
    print("\n# 4  near (1, 2), y1 - y2 = 1 + O(h^2)")
    exps = []
    for h in (1e-1, 1e-2, 1e-3, 1e-4):
        u, v = h, -h          # a generic direction, not the exact null direction
        y1 = (1 + u) ** 2 + (2 + v)
        y2 = (1 + u) * (2 + v)
        err = abs((y1 - y2) - 1.0)
        exps.append((h, err))
    # u^2 - u*v = 2h^2 for this direction: second order, constant 2.
    # tol 1e-6 absorbs the cancellation in (y1 - y2) - 1 at h = 1e-4.
    for h, err in exps:
        check("|y1-y2-1| / h^2 at h=%.0e" % h, err / (h * h), 2.0, 1e-6)
    # halving h must quarter the deviation
    r = exps[1][1] / exps[2][1]
    check("deviation ratio between h=1e-2 and 1e-3", r, 100.0, 2.0)


# ============================================ 5. two-link arm
def part5():
    print("\n# 5  planar two-link arm, l1 = l2 = 1")

    def arm_jac(t1, t2):
        s1, c1 = np.sin(t1), np.cos(t1)
        s12, c12 = np.sin(t1 + t2), np.cos(t1 + t2)
        return np.array([[-s1 - s12, -s12], [c1 + c12, c12]])

    t1, t2 = np.radians(18), np.radians(62)
    check("det J at 18/62 deg", np.linalg.det(arm_jac(t1, t2)), 0.8829, 5e-5)
    check("sin(theta2)", np.sin(t2), 0.8829, 5e-5)
    for t2 in (0.0, np.pi):
        check("det J at theta2=%.0f deg" % np.degrees(t2),
              np.linalg.det(arm_jac(t1, t2)), 0.0, 1e-12)
    # rank 1 at the singular pose: the velocity set is a segment
    J = arm_jac(np.radians(28), 0.0)
    check("rank at singular pose", np.linalg.matrix_rank(J), 1.0, 0)
    # the lost direction is along the link
    link = np.array([np.cos(np.radians(28)), np.sin(np.radians(28))])
    check("|J^T link| (link direction is lost)",
          np.linalg.norm(J.T @ link), 0.0, 1e-12)


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    print("\nALL CHECKS PASSED")
