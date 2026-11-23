#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Kronecker 符号与二阶张量运算的参考实现（教学用）。

范围：三维正交单位基底、实矢量与实二阶张量。矩阵与分量按列矢量约定：
A_ij 是 A 作用在第 j 个基矢量上的第 i 个分量，即 A_ij = e_i · (A e_j)。

记号：
- 重复指标求和（爱因斯坦约定），例如 a_i b_i = a·b。
- 张量积 (a⊗b)_ij = a_i b_j，对应矩阵外积 a bᵀ。
- 双点积 A:B = A_ij B_ij = tr(AᵀB)。
"""

import numpy as np


# ---------------------------------------------------------------- 基本分量
def delta(i, j):
    """Kronecker 符号 δ_ij（i、j 取 1..3）。"""
    return 1.0 if i == j else 0.0


def identity():
    """单位张量 I 的分量矩阵。"""
    return np.eye(3)


def tensor_product(a, b):
    """矢量张量积 a⊗b 的分量矩阵，即外积 a bᵀ。"""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return np.outer(a, b)


def basis_tensor(i, j):
    """基张量 e_i⊗e_j 的分量矩阵（只有 (i,j) 元为 1）。"""
    E = np.zeros((3, 3))
    E[i - 1, j - 1] = 1.0
    return E


def apply_tensor(A, v):
    """二阶张量作用于矢量：A v，分量 (Av)_i = A_ij v_j。"""
    return np.asarray(A, dtype=float) @ np.asarray(v, dtype=float)


def basis_action(i, j, v):
    """基张量作用：(e_i⊗e_j) v = v_j e_i。"""
    out = np.zeros(3)
    out[i - 1] = float(v[j - 1])
    return out


# ---------------------------------------------------------------- 乘法与转置
def multiply(A, B):
    """张量乘法（线性映射复合）C = A B 的分量矩阵，C_ij = A_ik B_kj。"""
    return np.asarray(A, dtype=float) @ np.asarray(B, dtype=float)


def transpose(A):
    """转置 Aᵀ，分量 (Aᵀ)_ij = A_ji。"""
    return np.asarray(A, dtype=float).T


# ---------------------------------------------------------------- 对称分解
def symmetric_part(A):
    """对称部分 S = (A + Aᵀ)/2。"""
    A = np.asarray(A, dtype=float)
    return 0.5 * (A + A.T)


def antisymmetric_part(A):
    """反对称部分 W = (A − Aᵀ)/2。"""
    A = np.asarray(A, dtype=float)
    return 0.5 * (A - A.T)


# ---------------------------------------------------------------- 迹与双点积
def trace(A):
    """迹 tr A = A_ii。"""
    return float(np.trace(np.asarray(A, dtype=float)))


def double_dot(A, B):
    """双点积 A:B = A_ij B_ij = tr(AᵀB)。"""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    return float(np.sum(A * B))


def double_dot_swapped(A, B):
    """交换指标的双点积 A_ij B_ji = tr(AB)。"""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    return float(np.sum(A * B.T))


def frobenius(A):
    """Frobenius 范数 ‖A‖_F = sqrt(A:A)。"""
    return float(np.sqrt(double_dot(A, A)))


# ---------------------------------------------------------------- 坐标变换
def rotation_matrix(axis, angle_rad):
    """绕单位轴 axis 转动 angle_rad 的右手正交矩阵。"""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    C = 1.0 - c
    return np.array([
        [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def rotate_tensor(A, Q):
    """正交基变换下的二阶张量：A' = Q A Qᵀ。"""
    A = np.asarray(A, dtype=float)
    Q = np.asarray(Q, dtype=float)
    return Q @ A @ Q.T


# ---------------------------------------------------------------- 算例数据
def article_vectors():
    """原文练习 2 的三个矢量。"""
    a = np.array([1.0, 2.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    c = np.array([3.0, 4.0, 5.0])
    return a, b, c


def article_decomposition():
    """原文 §7.3 的分解例。"""
    A = np.array([[1.0, 4.0, 0.0],
                  [2.0, 3.0, 0.0],
                  [0.0, 0.0, 5.0]])
    return A


def article_exercise3():
    """原文练习 3 的张量。"""
    return np.array([[1.0, 2.0, 0.0],
                     [0.0, 3.0, 0.0],
                     [0.0, 0.0, 4.0]])


def article_exercise4():
    """原文练习 4 的张量。"""
    return np.array([[2.0, 5.0, 0.0],
                     [1.0, 3.0, 0.0],
                     [0.0, 0.0, 4.0]])