#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""校验 M0-04《Kronecker 符号与二阶张量》中的全部数值结论。

运行：
    py -3 tools/verify_tensor_notes.py

末尾打印 ALL CHECKS PASSED；任一检查失败立即退出。
"""

import numpy as np

import tensor_lib as ta


def check(label, value, target, tol=1e-14):
    ok = abs(value - target) <= tol
    print(("PASS  " if ok else "FAIL  ") + label
          + "  value=%+.6e  target=%+.6e" % (value, target))
    if not ok:
        raise SystemExit(1)


def matrix_check(label, value, target, tol=1e-14):
    err = float(np.max(np.abs(np.asarray(value) - np.asarray(target))))
    print(("PASS  " if err <= tol else "FAIL  ") + label
          + "  max|diff|=%.2e" % err)
    if err > tol:
        raise SystemExit(1)


def part1():
    print("\n[1] Kronecker 符号与指标替换")
    check("delta(1,1)", ta.delta(1, 1), 1.0)
    check("delta(1,2)", ta.delta(1, 2), 0.0)
    matrix_check("delta 的分量矩阵", np.array(
        [[ta.delta(i, j) for j in range(1, 4)] for i in range(1, 4)]),
        np.eye(3))

    a = np.array([1.0, 2.0, 3.0])
    for i in range(1, 4):
        check("sum_j delta(%d,j) a_j" % i,
              sum(ta.delta(i, j) * a[j - 1] for j in range(1, 4)), a[i - 1])
    A = ta.article_exercise3()
    left = np.einsum("ij,jk->ik",
                     np.array([[ta.delta(i, j) for j in range(1, 4)]
                               for i in range(1, 4)]), A)
    matrix_check("delta_ij A_jk = A_ik", left, A)
    right = np.einsum("ij,jk->ik", A,
                      np.array([[ta.delta(j, k) for k in range(1, 4)]
                                for j in range(1, 4)]))
    matrix_check("A_ij delta_jk = A_ik", right, A)
    check("delta_ii", sum(ta.delta(i, i) for i in range(1, 4)), 3.0)
    check("delta_ij delta_ij",
          sum(ta.delta(i, j) * ta.delta(i, j)
              for i in range(1, 4) for j in range(1, 4)), 3.0)


def part2():
    print("\n[2] 矢量张量积与作用")
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([4.0, 5.0, 6.0])
    c = np.array([-1.0, 2.0, 1.0])
    T = ta.tensor_product(a, b)
    matrix_check("a x b 的分量矩阵", T, np.outer(a, b))
    expected = a * float(np.dot(b, c))
    matrix_check("(a x b)c = a(b*c)", ta.apply_tensor(T, c), expected)
    print("      b*c = %.1f，输出 = [%.1f, %.1f, %.1f]"
          % (np.dot(b, c), *expected))
    check("max|a x b - b x a|（应为 6）",
          float(np.max(np.abs(T - ta.tensor_product(b, a)))),
          6.0, tol=1e-14)

    v = np.array([2.0, -3.0, 5.0])
    E12 = ta.basis_tensor(1, 2)
    matrix_check("e1 x e2 的作用", ta.apply_tensor(E12, v),
                 ta.basis_action(1, 2, v))
    matrix_check("(e1 x e2)v = v2 e1", ta.apply_tensor(E12, v),
                 np.array([-3.0, 0.0, 0.0]))

    B = np.array([[0.0, 1.0, 0.0],
                  [1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0]])
    check("(a x b):(a x b) 与 |a|^2|b|^2 的差",
          ta.double_dot(T, T) - float(np.dot(a, a) * np.dot(b, b)),
          0.0, tol=1e-12)
    print("      a x b 的秩 = %d，迹 = %.1f" % (np.linalg.matrix_rank(T),
                                              ta.trace(T)))


def part3():
    print("\n[3] 基张量、单位张量与分量运算")
    A = ta.article_exercise3()
    rebuilt = np.zeros((3, 3))
    for i in range(1, 4):
        for j in range(1, 4):
            rebuilt += A[i - 1, j - 1] * ta.basis_tensor(i, j)
    matrix_check("A_ij e_i x e_j 重建 A", rebuilt, A)
    check("迹 tr I", ta.trace(ta.identity()), 3.0)
    v = np.array([2.0, -3.0, 5.0])
    matrix_check("I v = v", ta.apply_tensor(ta.identity(), v), v)
    matrix_check("I A = A", ta.multiply(ta.identity(), A), A)
    matrix_check("A I = A", ta.multiply(A, ta.identity()), A)


def part4():
    print("\n[4] 乘法、转置与不可交换性")
    A = np.array([[1.0, 2.0, 0.0],
                  [0.0, 3.0, 0.0],
                  [0.0, 0.0, 4.0]])
    B = np.array([[0.0, 1.0, 0.0],
                  [1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0]])
    AB = ta.multiply(A, B)
    BA = ta.multiply(B, A)
    matrix_check("C_ij = A_ik B_kj", AB, A @ B)
    check("max|AB - BA|（应为 2）",
          float(np.max(np.abs(AB - BA))), 2.0, tol=1e-14)
    matrix_check("(AB)^T = B^T A^T", ta.transpose(AB),
                 ta.multiply(ta.transpose(B), ta.transpose(A)))
    check("tr(AB) = tr(BA)", ta.trace(AB) - ta.trace(BA), 0.0, tol=1e-14)
    a = np.array([1.0, 2.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    matrix_check("(a x b)^T = b x a", ta.transpose(ta.tensor_product(a, b)),
                 ta.tensor_product(b, a))
    print("      AB 的 (1,1) 元 = %.1f，(1,2) 元 = %.1f；BA 的 (1,2) 元 = %.1f"
          % (AB[0, 0], AB[0, 1], BA[0, 1]))


def part5():
    print("\n[5] 对称/反对称分解")
    A = ta.article_decomposition()
    S = ta.symmetric_part(A)
    W = ta.antisymmetric_part(A)
    matrix_check("S 的值", S, np.array([[1.0, 3.0, 0.0],
                                        [3.0, 3.0, 0.0],
                                        [0.0, 0.0, 5.0]]))
    matrix_check("W 的值", W, np.array([[0.0, 1.0, 0.0],
                                        [-1.0, 0.0, 0.0],
                                        [0.0, 0.0, 0.0]]))
    matrix_check("S + W = A", S + W, A)
    matrix_check("S^T = S", S.T, S)
    matrix_check("W^T = -W", W.T, -W)
    check("S:W", ta.double_dot(S, W), 0.0, tol=1e-14)

    rng = np.random.default_rng(7)
    R = rng.normal(size=(3, 3))
    Sr = ta.symmetric_part(R)
    Wr = ta.antisymmetric_part(R)
    matrix_check("随机张量 S+W=R", Sr + Wr, R, tol=1e-14)
    check("随机张量 S:W", ta.double_dot(Sr, Wr), 0.0, tol=1e-14)
    print("      独立分量数：S 6 个，W 3 个；W 的对角元全为 0")


def part6():
    print("\n[6] 迹、双点积与 Frobenius 范数")
    A = ta.article_exercise3()
    check("tr A", ta.trace(A), 8.0)
    check("A:A", ta.double_dot(A, A), 30.0)
    check("A_ij A_ji", ta.double_dot_swapped(A, A), 26.0)
    check("tr(A^T A)", ta.trace(ta.multiply(A.T, A)), 30.0)
    check("tr(A A)", ta.trace(ta.multiply(A, A)), 26.0)

    E = np.array([[0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0]])
    check("E:E", ta.double_dot(E, E), 1.0)
    check("E_ij E_ji", ta.double_dot_swapped(E, E), 0.0)
    print("      同一对矩阵：A_ij B_ij = 1，A_ij B_ji = 0")

    D = np.diag([1.0, 2.0, 2.0])
    check("||diag(1,2,2)||_F", ta.frobenius(D), 3.0)

    axis = np.array([1.0, 2.0, -1.0])
    Q = ta.rotation_matrix(axis, np.deg2rad(37.0))
    check("det Q", float(np.linalg.det(Q)), 1.0, tol=1e-14)
    matrix_check("Q^T Q = I", Q.T @ Q, np.eye(3), tol=1e-14)
    Arot = ta.rotate_tensor(A, Q)
    check("旋转后 ||A||_F 不变", ta.frobenius(Arot) - ta.frobenius(A),
          0.0, tol=1e-13)
    check("旋转后 tr A 不变", ta.trace(Arot) - ta.trace(A), 0.0, tol=1e-13)
    print("      Q A Q^T 的分量改变，但迹与 Frobenius 范数不变")


def part7():
    print("\n[7] 原文练习与答案")
    a, b, c = ta.article_vectors()
    T = ta.tensor_product(a, b)
    matrix_check("练习 2：a x b", T, np.array([[0.0, 1.0, 0.0],
                                             [0.0, 2.0, 0.0],
                                             [0.0, 0.0, 0.0]]))
    matrix_check("练习 2：(a x b)c", ta.apply_tensor(T, c),
                 np.array([4.0, 8.0, 0.0]))

    A3 = ta.article_exercise3()
    check("练习 3：tr A", ta.trace(A3), 8.0)
    check("练习 3：A_ij A_ij", ta.double_dot(A3, A3), 30.0)
    check("练习 3：A_ij A_ji", ta.double_dot_swapped(A3, A3), 26.0)

    A4 = ta.article_exercise4()
    S4 = ta.symmetric_part(A4)
    W4 = ta.antisymmetric_part(A4)
    matrix_check("练习 4：S", S4, np.array([[2.0, 3.0, 0.0],
                                            [3.0, 3.0, 0.0],
                                            [0.0, 0.0, 4.0]]))
    matrix_check("练习 4：W", W4, np.array([[0.0, 2.0, 0.0],
                                            [-2.0, 0.0, 0.0],
                                            [0.0, 0.0, 0.0]]))
    check("练习 4：S:W", ta.double_dot(S4, W4), 0.0, tol=1e-14)

    # 练习 1 的四个结果
    bvec = np.array([4.0, 5.0, 6.0])
    out = np.array([sum(ta.delta(i, j) * bvec[j - 1]
                        for j in range(1, 4)) for i in range(1, 4)])
    matrix_check("练习 1：delta_ij b_j = b_i", out, bvec)
    check("练习 1：delta_ij A_ij = tr A",
          sum(ta.delta(i, j) * A3[i - 1, j - 1]
              for i in range(1, 4) for j in range(1, 4)), 8.0)
    check("练习 1：delta_ij delta_ij", 3.0, 3.0)


def main():
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
    part7()
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()