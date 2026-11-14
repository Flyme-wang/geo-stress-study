#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用户子程序 FRIC 的线性黏性摩擦参考实现（教学用）。

Abaqus 的 FRIC 允许用户用子程序覆盖库仑摩擦定律。子程序在每个接触点、
每个摩擦方向上被调用，输入增量滑移 dgam、总滑移 slip、两侧温度 temp、
时间增量 dtime 等，返回剪应力 tau 及其切线。

本模块只实现原文两个算例所用的线性黏性定律

    tau = dgam / (k * dtime)       等价于       dgam = k * tau * dtime

Abaqus 的 dgam / slip 是**长度**，所以严格写应是滑移速度

    hdot = k * tau

其中 k 的单位是“长度/(应力*时间)”，例如 m/(Pa*s) 或 mm/(Pa*s)；
长度的单位跟模型一致，k 的单位也随之改变。原文写 gamma_dot = k*tau
并称 gamma_dot 为剪应变率，但 Abaqus 的 dgam 不是应变，这是需要纠正的口径。
若按单位面积计，则“力”和“应力”的数值相同；面积不是 1 时必须先除以面积。
"""

import numpy as np


def slip_increment(tau, k, dt):
    """线性黏性界面的增量滑移 dgam = k*tau*dt（长度单位随模型）。"""
    return k * tau * dt


def slip_rate(tau, k):
    """线性黏性界面的滑移速度 dgam/dt = k*tau（长度/时间）。"""
    return k * tau


# ---------------------------------------------------------------- 算例一
def example1_steps(force_mag=100.0, k=0.001, dt=1.0):
    """B31 梁压在刚性平面上，三个增量依次沿 x、45°、y 施加切向力。

    返回每个增量的力矢量、剪应力矢量、增量滑移矢量与累计滑移矢量。
    按单位面积计，因此力矢量与剪应力矢量数值相同。
    """
    directions = [
        ("沿 x", np.array([1.0, 0.0])),
        ("45°", np.array([1.0, 1.0]) / np.sqrt(2.0)),
        ("沿 y", np.array([0.0, 1.0])),
    ]
    cum = np.zeros(2)
    path = 0.0
    rows = []
    for name, unit in directions:
        force = force_mag * unit
        dgam = slip_increment(force, k, dt)
        cum = cum + dgam
        path += float(np.linalg.norm(dgam))
        rows.append({
            "name": name,
            "force": force,
            "tau": force.copy(),
            "dgam": dgam,
            "cum": cum.copy(),
            "path": path,
        })
    return rows


def example1_summary(force_mag=100.0, k=0.001, dt=1.0):
    """算例一的终点、路径长、弦长与方向角。"""
    rows = example1_steps(force_mag, k, dt)
    final = rows[-1]["cum"]
    path = rows[-1]["path"]
    chord = float(np.linalg.norm(final))
    angle = float(np.degrees(np.arctan2(final[1], final[0])))
    return {
        "final": final,
        "path": path,
        "chord": chord,
        "angle": angle,
    }


# ---------------------------------------------------------------- 算例二
def example2_temperature(t):
    """B 块温度从 0 线性升到 200，A 块保持 0；返回 (theta_A, theta_B)。"""
    return 0.0, 200.0 * t


def example2_k(t, k0=0.001, k1=1.0e-5):
    """k = k0 + k1*(theta_A+theta_B)/2 = k0 + k1*100*t。"""
    theta_a, theta_b = example2_temperature(t)
    return k0 + k1 * 0.5 * (theta_a + theta_b)


def example2_continuous(k0=0.001, k1=1.0e-5, tau=100.0):
    """连续解：gamma = tau * integral_0^1 k(t) dt = 0.15。"""
    # k(t) = k0 + k1*100*t，积分 = k0 + 0.5*k1*100
    return tau * (k0 + 0.5 * k1 * 100.0)


def example2_end_slip(n, k0=0.001, k1=1.0e-5, tau=100.0):
    """每步用增量末温度采样时，n 步后的总滑移。"""
    i = np.arange(1, n + 1, dtype=float)
    t = i / n
    k = k0 + k1 * 100.0 * t
    return float(np.sum(k * tau * (1.0 / n)))


def example2_mid_slip(n, k0=0.001, k1=1.0e-5, tau=100.0):
    """每步用增量中点温度采样时，n 步后的总滑移（对线性 k 精确）。"""
    i = np.arange(n, dtype=float)
    t = (i + 0.5) / n
    k = k0 + k1 * 100.0 * t
    return float(np.sum(k * tau * (1.0 / n)))


def example2_end_slip_closed(n, k0=0.001, k1=1.0e-5, tau=100.0):
    """增量末采样的闭式：0.15 + 0.05/n（本文参数下）。"""
    return example2_continuous(k0, k1, tau) + 0.05 / n


def example2_mid_slip_closed(n, k0=0.001, k1=1.0e-5, tau=100.0):
    """增量中点采样的闭式：与 n 无关，恒为 0.15。"""
    return example2_continuous(k0, k1, tau)


# ---------------------------------------------------------------- 切线
def fric_tangents(dgam, dt, k0=0.001, k1=1.0e-5, theta_bar=0.0, k=None):
    """返回 (tau, dtau/ddgam, dtau/dtheta1=dtau/dtheta2)。

    由 tau = dgam/(k*dt)，k = k0 + k1*(theta1+theta2)/2：
        dtau/ddgam = 1/(k*dt)
        dtau/dtheta1 = -k1*dgam/(2*k^2*dt)
    """
    if k is None:
        k = k0 + k1 * theta_bar
    tau = dgam / (k * dt)
    ddtddg = 1.0 / (k * dt)
    ddtddt = -k1 * dgam / (2.0 * k * k * dt)
    return tau, ddtddg, ddtddt


# ---------------------------------------------------------------- 摩擦热
def friction_heat(k, tau, eta=0.5):
    """单位面积摩擦热流 q = eta*k*tau^2（W/m^2）。"""
    return eta * k * tau * tau


def example2_heat_continuous(k0=0.001, k1=1.0e-5, tau=100.0, eta=0.5):
    """连续解的累计摩擦热（J/m^2）：0.5*100^2*0.0015 = 7.5。"""
    integral_k = k0 + 0.5 * k1 * 100.0
    return eta * tau * tau * integral_k


def example2_heat_end(n, k0=0.001, k1=1.0e-5, tau=100.0, eta=0.5):
    """增量末采样、n 步累计的摩擦热。"""
    i = np.arange(1, n + 1, dtype=float)
    t = i / n
    k = k0 + k1 * 100.0 * t
    return float(np.sum(eta * k * tau * tau * (1.0 / n)))


def example2_heat_mid(n, k0=0.001, k1=1.0e-5, tau=100.0, eta=0.5):
    """增量中点采样、n 步累计的摩擦热（对线性 k 精确）。"""
    i = np.arange(n, dtype=float)
    t = (i + 0.5) / n
    k = k0 + k1 * 100.0 * t
    return float(np.sum(eta * k * tau * tau * (1.0 / n)))


# ---------------------------------------------------------------- 库仑对照
def coulomb_limit(mu, p):
    """库仑摩擦上限 tau_max = mu*p。"""
    return mu * p


def coulomb_regime(tau, mu, p):
    """返回 'stick' 或 'slip'（临界点按 slip 处理）。"""
    return "stick" if abs(tau) < coulomb_limit(mu, p) else "slip"
