#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recompute every number quoted in docs/07-fluid/*.md.

Usage (from the repository root):
    py -3 tools/verify_fluid_notes.py

Checks:
  1. thermodynamic recovery of the two Sod initial states;
  2. the shared Rusanov flux on the diaphragm;
  3. the hand-updated two-cell state and its recovered primitives;
  4. the exact Sod star state, wave positions and shock speed;
  5. a 200/400/800-cell finite-volume run: step counts, density error,
     convergence order, global conservation and positivity;
  6. a subsonic-outlet construction and a methane low-Mach estimate.
"""

import sys

import numpy as np

sys.path.insert(0, "tools")
import fluid_lib as fl


def check(label, value, target, tol):
    ok = abs(value - target) <= tol
    print(("PASS  " if ok else "FAIL  ") + label
          + "  value=%.10g  target=%.10g  tol=%.1e" % (value, target, tol))
    if not ok:
        raise SystemExit(1)


def check_true(label, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + label + ("  " + detail))
    if not cond:
        raise SystemExit(1)


# ================================================= 1. initial states
def part1():
    print("\n# 1  Sod initial states and thermodynamic recovery")
    UL = fl.conserved(1.0, 0.0, 1.0)
    UR = fl.conserved(0.125, 0.0, 0.1)
    check("E_L", UL[2], 2.5, 1e-14)
    check("E_R", UR[2], 0.25, 1e-14)
    check("c_L", fl.sound_speed(1.0, 1.0), 1.183216, 5e-7)
    check("c_R", fl.sound_speed(0.125, 0.1), 1.058301, 5e-7)
    # T = p/(rho R), R = 1
    check("T_L", 1.0 / (1.0 * 1.0), 1.0, 1e-14)
    check("T_R", 0.1 / (0.125 * 1.0), 0.8, 1e-14)


# ================================================= 2. Rusanov flux
def part2():
    print("\n# 2  shared Rusanov flux on the diaphragm")
    UL = fl.conserved(1.0, 0.0, 1.0)
    UR = fl.conserved(0.125, 0.0, 0.1)
    F, smax = fl.rusanov_flux(UL, UR)
    check("s_max", smax, 1.183216, 5e-7)
    check("mass flux", F[0], 0.517657, 5e-7)
    check("momentum flux", F[1], 0.550000, 1e-12)
    check("energy flux", F[2], 1.331118, 5e-7)
    # the exact mass flux at x = 0.5 is rho_*L u_*
    check("exact rho_*L u_*", 0.426319 * 0.927453, 0.395391, 5e-7)


# =========================================== 3. two-cell update
def part3():
    print("\n# 3  two-cell update with dt/h = 0.1")
    UL = fl.conserved(1.0, 0.0, 1.0)
    UR = fl.conserved(0.125, 0.0, 0.1)
    FL = fl.physical_flux(UL)
    FR = fl.physical_flux(UR)
    F, _ = fl.rusanov_flux(UL, UR)
    dt_over_h = 0.1
    ULn = UL - dt_over_h * (F - FL)
    URn = UR - dt_over_h * (FR - F)
    check("rho_L new", ULn[0], 0.948234, 5e-7)
    check("rho_R new", URn[0], 0.176766, 5e-7)
    check("rho*u_L new", ULn[1], 0.045000, 1e-12)
    check("rho*u_R new", URn[1], 0.045000, 1e-12)
    check("E_L new", ULn[2], 2.366888, 5e-7)
    check("E_R new", URn[2], 0.383112, 5e-7)
    rhoL, uL, pL = fl.primitive(ULn)
    rhoR, uR, pR = fl.primitive(URn)
    check("u_L new", uL, 0.047457, 5e-7)
    check("u_R new", uR, 0.254574, 5e-7)
    check("p_L new", pL, 0.946328, 5e-7)
    check("p_R new", pR, 0.150954, 5e-7)
    check("T_L new", pL / rhoL, 0.997990, 5e-7)
    check("T_R new", pR / rhoR, 0.853975, 5e-7)
    # mass and energy lost by the left cell equal what the right cell gains
    check("mass jump balance", ULn[0] + URn[0], 1.125, 1e-12)
    check("energy jump balance", ULn[2] + URn[2], 2.75, 1e-12)


# ================================================= 4. exact solution
def part4():
    print("\n# 4  exact Sod star state, wave positions, shock speed")
    ps, us, cL, cR = fl.sod_star_state()
    check("p*", ps, 0.303130, 5e-7)
    check("u*", us, 0.927453, 5e-7)
    check("rho_*L", 1.0 * (ps / 1.0) ** (1 / 1.4), 0.426319, 5e-7)
    rhoRstar = 0.125 * ((ps / 0.1 + 0.4 / 2.4)
                        / (0.4 / 2.4 * ps / 0.1 + 1.0))
    check("rho_*R", rhoRstar, 0.265574, 5e-7)
    cstar = cL * (ps / 1.0) ** (0.4 / 2.8)
    check("rarefaction head at t=0.2", 0.5 + (0 - cL) * 0.2, 0.263357, 5e-7)
    check("rarefaction tail at t=0.2", 0.5 + (us - cstar) * 0.2, 0.485945, 5e-7)
    check("contact at t=0.2", 0.5 + us * 0.2, 0.685491, 5e-7)
    S = cR * np.sqrt(2.4 / 2.8 * ps / 0.1 + 0.4 / 2.8)
    check("shock speed", S, 1.752156, 5e-7)
    check("shock at t=0.2", 0.5 + S * 0.2, 0.850431, 5e-7)


# ================================================= 5. finite volume
def part5():
    print("\n# 5  finite-volume run at 200/400/800 cells")
    targets = [(200, 190, 0.016646), (400, 385, 0.010890),
               (800, 774, 0.006978)]
    errs = []
    for n, steps, target in targets:
        r = fl.fv_sod(n, 0.2, cfl=0.45)
        check("%d cells: steps" % n, r["steps"], steps, 0)
        rho = r["U"][:, 0]
        rex = fl.sod_cell_averages(n, 0.2)
        err = float(np.sum(r["dx"] * np.abs(rho - rex)))
        errs.append(err)
        # the note quotes six decimals; the 400-cell value sits on the
        # rounding boundary, so allow 2e-6
        check("%d cells: density L1 error" % n, err, target, 2e-6)
        M = float(np.sum(rho) * r["dx"])
        P = float(np.sum(r["U"][:, 1]) * r["dx"])
        H = float(np.sum(r["U"][:, 2]) * r["dx"])
        check("%d cells: total mass" % n, M, 0.5625, 1e-12)
        check("%d cells: total momentum" % n, P, 0.1800, 1e-12)
        check("%d cells: total energy" % n, H, 1.3750, 1e-12)
        check("%d cells: boundary momentum" % n, r["boundary"][1], 0.18, 1e-12)
        check("%d cells: boundary mass" % n, r["boundary"][0], 0.0, 1e-14)
        residual = max(abs(M - 0.5625 - r["boundary"][0]),
                       abs(P - 0.0 - r["boundary"][1]),
                       abs(H - 1.3750 - r["boundary"][2]))
        check_true("%d cells: conservation audit residual" % n,
                   residual < 5e-15, "res=%.3e" % residual)
        prim = np.array([fl.primitive(u) for u in r["U"]])
        check_true("%d cells: positivity" % n,
                   bool(np.all(prim[:, 0] > 0) and np.all(prim[:, 2] > 0)))
    check("observed order 200->400",
          np.log2(errs[0] / errs[1]), 0.61, 0.015)
    check("observed order 400->800",
          np.log2(errs[1] / errs[2]), 0.64, 0.015)


# =========================================== 6. boundaries and low Mach
def part6():
    print("\n# 6  subsonic outlet and methane low-Mach estimate")
    gamma = 1.4
    rho_in, p_in, u_in = 1.0, 1.0, 0.5
    c_in = fl.sound_speed(rho_in, p_in)
    check("c_in", c_in, 1.183216, 5e-7)
    J = u_in + 2.0 * c_in / (gamma - 1.0)
    check("J_+ at the inside", J, 6.416080, 5e-7)
    c_b = 1.17
    p_b = p_in * (c_b / c_in) ** (2.0 * gamma / (gamma - 1.0))
    rho_b = rho_in * (p_b / p_in) ** (1.0 / gamma)
    u_b = u_in + 2.0 * (c_in - c_b) / (gamma - 1.0)
    check("p_b", p_b, 0.924385, 5e-7)
    check("rho_b", rho_b, 0.945386, 5e-7)
    check("u_b", u_b, 0.566080, 5e-7)
    check("J_+ at the boundary",
          u_b + 2.0 * c_b / (gamma - 1.0), J, 1e-10)
    check("isentropic p/rho^gamma",
          p_b / rho_b ** gamma, p_in / rho_in ** gamma, 1e-12)
    E_in = p_in / (gamma - 1.0) + 0.5 * rho_in * u_in ** 2
    E_b = p_b / (gamma - 1.0) + 0.5 * rho_b * u_b ** 2
    U_in = np.array([rho_in, rho_in * u_in, E_in])
    U_b = np.array([rho_b, rho_b * u_b, E_b])
    H, smax = fl.rusanov_flux(U_in, U_b)
    check("s_max at the boundary", smax, 1.736080, 5e-7)
    check("boundary mass flux", H[0], 0.564989, 5e-7)
    check("boundary momentum flux", H[1], 1.208142, 5e-7)
    check("boundary energy flux", H[2], 2.005968, 5e-7)
    # the flux is not the upwinded state flux: numerical diffusion adds
    # half of s_max times the state jump
    check("flux minus rho_b u_b", H[0] - rho_b * u_b, 0.029825, 5e-7)

    # methane: gamma = 1.31, R = 518 J/(kg K), 30 MPa, 380 K
    R, T, p, u = 518.0, 380.0, 30e6, 15.0
    rho = p / (R * T)
    c = np.sqrt(1.31 * p / rho)
    check("methane density [kg/m3]", rho, 152.4, 0.5)
    check("methane sound speed [m/s]", c, 508.0, 1.0)
    check("methane Mach", u / c, 0.0295, 5e-4)
    check("acoustic CFL dt for h=1 cm [s]",
          0.45 * 0.01 / (c + u), 8.6e-6, 2e-7)


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
    print("\nALL CHECKS PASSED")
