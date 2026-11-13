#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""一维可压缩 Euler 的小型参考实现（教学用，不追求性能）。

提供三件事：
1. 定比热理想气体的守恒量 <-> 原始变量换算；
2. Rusanov（局部 Lax-Friedrichs）数值通量；
3. Sod 激波管的精确 Riemann 解与一维有限体积解。

所有量均为无量纲：R = 1, gamma = 1.4, 计算域 0 <= x <= 1。
"""

import numpy as np
from scipy.optimize import brentq

GAMMA = 1.4


def primitive(U, gamma=GAMMA):
    """守恒量 U=(rho, rho*u, E) -> (rho, u, p)。"""
    rho = U[0]
    m = U[1]
    E = U[2]
    u = m / rho
    p = (gamma - 1.0) * (E - 0.5 * m * m / rho)
    return rho, u, p


def conserved(rho, u, p, gamma=GAMMA):
    """(rho, u, p) -> 守恒量 U=(rho, rho*u, E)。"""
    E = p / (gamma - 1.0) + 0.5 * rho * u * u
    return np.array([rho, rho * u, E], dtype=float)


def sound_speed(rho, p, gamma=GAMMA):
    return np.sqrt(gamma * p / rho)


def physical_flux(U, gamma=GAMMA):
    """一维 Euler 物理通量 F=(rho*u, rho*u^2+p, (E+p)*u)。"""
    rho, u, p = primitive(U, gamma)
    E = U[2]
    return np.array([rho * u, rho * u * u + p, (E + p) * u], dtype=float)


def rusanov_flux(UL, UR, gamma=GAMMA):
    """Rusanov 共享通量，返回 (Fhat, smax)。"""
    FL = physical_flux(UL, gamma)
    FR = physical_flux(UR, gamma)
    rhoL, uL, pL = primitive(UL, gamma)
    rhoR, uR, pR = primitive(UR, gamma)
    smax = max(abs(uL) + sound_speed(rhoL, pL, gamma),
               abs(uR) + sound_speed(rhoR, pR, gamma))
    return 0.5 * (FL + FR) - 0.5 * smax * (UR - UL), smax


# ------------------------------------------------------------ Sod 精确解
def _f_k(p, pK, rhoK, cK, gamma):
    """压力函数：激波分支（p>pK）或稀疏分支（p<=pK）。"""
    if p > pK:
        A = 2.0 / ((gamma + 1.0) * rhoK)
        B = (gamma - 1.0) / (gamma + 1.0) * pK
        return (p - pK) * np.sqrt(A / (p + B))
    return 2.0 * cK / (gamma - 1.0) * (
        (p / pK) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)


def sod_star_state(gamma=GAMMA, rhoL=1.0, uL=0.0, pL=1.0,
                   rhoR=0.125, uR=0.0, pR=0.1):
    """求中间压力 p* 与中间速度 u*。"""
    cL = sound_speed(rhoL, pL, gamma)
    cR = sound_speed(rhoR, pR, gamma)

    def f(p):
        return (_f_k(p, pL, rhoL, cL, gamma)
                + _f_k(p, pR, rhoR, cR, gamma) + uR - uL)

    lo = 1e-12
    hi = max(pL, pR)
    while f(hi) < 0:
        hi *= 2.0
    ps = brentq(f, lo, hi, xtol=1e-14, rtol=1e-14)
    us = 0.5 * (uL + uR) + 0.5 * (
        _f_k(ps, pR, rhoR, cR, gamma) - _f_k(ps, pL, rhoL, cL, gamma))
    return ps, us, cL, cR


def sod_exact(x, t, gamma=GAMMA, rhoL=1.0, uL=0.0, pL=1.0,
              rhoR=0.125, uR=0.0, pR=0.1, x0=0.5):
    """Sod 激波管精确解在 (x,t) 处的原始变量。

    支持标量或数组 x；t 必须为正。
    """
    ps, us, cL, cR = sod_star_state(gamma, rhoL, uL, pL, rhoR, uR, pR)
    x = np.atleast_1d(np.asarray(x, dtype=float))
    xi = (x - x0) / t

    rho = np.empty_like(xi)
    u = np.empty_like(xi)
    p = np.empty_like(xi)

    # 左波：Sod 中 p* < pL，是稀疏波
    if ps < pL:
        cstarL = cL * (ps / pL) ** ((gamma - 1.0) / (2.0 * gamma))
        rhoLstar = rhoL * (ps / pL) ** (1.0 / gamma)
        head = uL - cL
        tail = us - cstarL
        fan = (xi >= head) & (xi < tail)
        u[fan] = 2.0 / (gamma + 1.0) * (cL + 0.5 * (gamma - 1.0) * uL + xi[fan])
        cfan = 2.0 / (gamma + 1.0) * (
            cL + 0.5 * (gamma - 1.0) * (uL - xi[fan]))
        rho[fan] = rhoL * (cfan / cL) ** (2.0 / (gamma - 1.0))
        p[fan] = pL * (cfan / cL) ** (2.0 * gamma / (gamma - 1.0))
    else:
        rhoLstar = rhoL * (
            (ps / pL + (gamma - 1.0) / (gamma + 1.0))
            / ((gamma - 1.0) / (gamma + 1.0) * ps / pL + 1.0))
        tail = us

    # 右波：Sod 中 p* > pR，是激波
    if ps > pR:
        rhoRstar = rhoR * (
            (ps / pR + (gamma - 1.0) / (gamma + 1.0))
            / ((gamma - 1.0) / (gamma + 1.0) * ps / pR + 1.0))
        S = uR + cR * np.sqrt(
            (gamma + 1.0) / (2.0 * gamma) * ps / pR
            + (gamma - 1.0) / (2.0 * gamma))
    else:
        rhoRstar = rhoR * (ps / pR) ** (1.0 / gamma)
        S = uR + cR

    left = xi < us
    # 左常值区
    m = left & (xi < (uL - cL))
    rho[m], u[m], p[m] = rhoL, uL, pL
    # 左星区
    m = left & (xi >= tail)
    rho[m], u[m], p[m] = rhoLstar, us, ps
    # 右星区
    m = (~left) & (xi <= S)
    rho[m], u[m], p[m] = rhoRstar, us, ps
    # 右常值区
    m = (~left) & (xi > S)
    rho[m], u[m], p[m] = rhoR, uR, pR
    return rho, u, p


def sod_cell_averages(n, t, gamma=GAMMA, nsub=401):
    """精确解在 n 个等宽单元上的平均值（用于收敛阶）。"""
    dx = 1.0 / n
    xs = (np.arange(n)[:, None] + (np.arange(nsub)[None, :] + 0.5) / nsub) * dx
    rho, u, p = sod_exact(xs.ravel(), t, gamma)
    return rho.reshape(n, nsub).mean(axis=1)


def fv_sod(n, t_end, cfl=0.45, gamma=GAMMA):
    """一维有限体积：分片常数 + Rusanov + 前向 Euler，外推边界。

    返回 dict：x（单元中心）、U（末态）、步数、边界收支累计。
    """
    dx = 1.0 / n
    x = (np.arange(n) + 0.5) * dx
    U = np.empty((n, 3))
    UL = conserved(1.0, 0.0, 1.0, gamma)
    UR = conserved(0.125, 0.0, 0.1, gamma)
    U[x < 0.5] = UL
    U[x >= 0.5] = UR

    t = 0.0
    steps = 0
    inflow_mass = 0.0
    inflow_mom = 0.0
    inflow_energy = 0.0
    while t < t_end - 1e-15:
        smax = 0.0
        for i in range(n):
            rho, u, p = primitive(U[i], gamma)
            smax = max(smax, abs(u) + sound_speed(rho, p, gamma))
        dt = cfl * dx / smax
        if t + dt > t_end:
            dt = t_end - t
        # 面通量：内部共享；两端外推（ghost = 边值）
        F = np.empty((n + 1, 3))
        F[0] = physical_flux(U[0], gamma)
        F[-1] = physical_flux(U[-1], gamma)
        for i in range(1, n):
            F[i], _ = rusanov_flux(U[i - 1], U[i], gamma)
        Unew = U - dt / dx * (F[1:] - F[:-1])
        # 边界交换（外法向）：左端流入为正，右端流出为正
        inflow_mass += dt * (F[0, 0] - F[-1, 0])
        inflow_mom += dt * (F[0, 1] - F[-1, 1])
        inflow_energy += dt * (F[0, 2] - F[-1, 2])
        U = Unew
        t += dt
        steps += 1
    return {
        "x": x, "U": U, "steps": steps,
        "boundary": np.array([inflow_mass, inflow_mom, inflow_energy]),
        "dx": dx,
    }
