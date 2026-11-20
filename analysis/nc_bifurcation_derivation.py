"""Locate a Hopf point in the fast subsystem and test its criticality.

The slow variables are held fixed at each operating point. The first
Lyapunov coefficient distinguishes the two signs of criticality under the
implemented convention. A separate integration checks whether onset is
consistent with the supercritical amplitude law. Both are numerical results
for the specified model, with no experimental data involved.
"""
from __future__ import annotations
import json, os
import numpy as np
from scipy.optimize import brentq

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "model"))
from params import P  # frozen v0 params

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "nc")
A, B_, EPS = P["a_fhn"], P["b_fhn"], P["eps_fhn"]


def _kappa(c, E):
    return (E / (P["K_E"] + E)) / P["tau_fast"]


def _gstar(tau):
    """Fixed point G*: dw=0 -> w=(G+a)/b; dG=0 -> G - G^3/3 - (G+a)/b + tau = 0."""
    f = lambda G: G - G**3 / 3.0 - (G + A) / B_ + tau
    # bracket a root over the FHN range
    xs = np.linspace(-2.5, 2.5, 400)
    fv = [f(x) for x in xs]
    roots = []
    for i in range(len(xs) - 1):
        if fv[i] == 0 or fv[i] * fv[i + 1] < 0:
            roots.append(brentq(f, xs[i], xs[i + 1]))
    return roots


def _jac(G, kappa):
    return np.array([[kappa * (1 - G**2), -kappa], [EPS, -EPS * B_]])


def _l1(G, kappa, omega0):
    """First Lyapunov coefficient (Kuznetsov 2D). Nonlinearity only in dG: F_G = -kappa*G^3/3.
    B(x,y)_1 = F_GG*x1*y1 = -2*kappa*G*x1*y1 ; C(x,y,z)_1 = F_GGG*x1*y1*z1 = -2*kappa*x1*y1*z1."""
    Aj = _jac(G, kappa)
    wv, V = np.linalg.eig(Aj)
    k = int(np.argmin(np.abs(wv.imag - omega0)))     # eigenvalue +i*omega0
    q = V[:, k]
    # adjoint eigenvector p: A^T p = -i omega0 p ; normalise <p,q>=1 (conj-linear first arg)
    wpv, Vp = np.linalg.eig(Aj.T)
    kp = int(np.argmin(np.abs(wpv - np.conj(wv[k]))))
    p = Vp[:, kp]
    q = q / np.sqrt(np.vdot(q, q).real)
    p = p / np.conj(np.vdot(p, q))                   # <p,q> = conj(p).q = 1
    def Bf(x, y): return np.array([-2 * kappa * G * x[0] * y[0], 0.0], complex)
    def Cf(x, y, z): return np.array([-2 * kappa * x[0] * y[0] * z[0], 0.0], complex)
    g20 = np.vdot(p, Bf(q, q))
    g11 = np.vdot(p, Bf(q, np.conj(q)))
    g21 = np.vdot(p, Cf(q, q, np.conj(q)))
    return float((1j * g20 * g11 + omega0 * g21).real / (2 * omega0**2))


def find_hopf(c, E, tau_grid=np.linspace(-1.5, 1.5, 3001)):
    kappa = _kappa(c, E)
    prev = None
    for tau in tau_grid:
        for G in _gstar(tau):
            J = _jac(G, kappa)
            ev = np.linalg.eigvals(J)
            tr = ev.real.max() if np.iscomplexobj(ev) else np.trace(J)
            # focus with complex pair; track the real part of the complex eigenvalue
            if abs(ev[0].imag) > 1e-9:
                re = ev[0].real
                if prev is not None and prev[1] * re < 0:   # sign change -> Hopf crossing
                    tau_c = 0.5 * (prev[0] + tau); G_c = _gstar(tau_c)
                    G_c = min(G_c, key=lambda g: abs(_jac(g, kappa)[0, 0] + _jac(g, kappa)[1, 1])) if G_c else G
                    omega0 = float(np.sqrt(max(np.linalg.det(_jac(G_c, kappa)), 1e-12)))
                    l1 = _l1(G_c, kappa, omega0)
                    return dict(tau_c=float(tau_c), G_c=float(G_c), omega0=omega0,
                                period_norm=float(2 * np.pi / omega0), l1=l1,
                                supercritical=bool(l1 < 0), kappa=float(kappa))
                prev = (tau, re)
    return None


def onset_scaling(hopf, c, E, n=9):
    """Integrate the fast subsystem just ABOVE tau_c; limit-cycle amplitude should scale ~sqrt(tau-tau_c)."""
    from scipy.integrate import solve_ivp
    kappa = _kappa(c, E); tc = hopf["tau_c"]
    def rhs(t, x, tau):
        G, w = x
        return [kappa * (G - G**3 / 3 - w + tau), EPS * (G + A - B_ * w)]
    dtaus = np.linspace(0.01, 0.20, n); amps = []
    for d in dtaus:
        tau = tc + d
        gs = _gstar(tau); G0 = (gs[0] if gs else 0.0) + 0.05
        sol = solve_ivp(rhs, (0, 4000), [G0, (G0 + A) / B_], args=(tau,), rtol=1e-8, atol=1e-10, max_step=2.0)
        g = sol.y[0][sol.t > 3000]
        amps.append(float((g.max() - g.min()) / 2))
    amps = np.array(amps)
    # supercritical: amp^2 linear in (tau-tau_c) through ~origin -> fit amp^2 = m*dtau + k
    m, k = np.polyfit(dtaus, amps**2, 1)
    r2 = 1 - np.sum((amps**2 - (m * dtaus + k))**2) / (np.sum((amps**2 - np.mean(amps**2))**2) + 1e-12)
    return dict(dtaus=dtaus.tolist(), amps=amps.tolist(), sqrt_law_slope=float(m),
                sqrt_law_intercept=float(k), sqrt_law_r2=float(r2),
                consistent_supercritical=bool(m > 0 and r2 > 0.9 and abs(k) < 0.05))


def main():
    os.makedirs(OUT, exist_ok=True)
    c0, E0 = P["c0"], P["E0"]
    base = find_hopf(c0, E0)
    res = {"provenance": "analysis-of-frozen-model (D2); frozen v0 params; deterministic",
           "control_parameter": "RhoA-driven contractile tone tau = I + k_c*(1-c)",
           "operating_point": {"c": c0, "E": E0}, "hopf": base}
    if base:
        res["onset"] = onset_scaling(base, c0, E0)
    # At fixed tone, c enters only through its contribution to that tone.
    # Vary energy alone to check distinct operating points.
    grid = []
    for E in np.linspace(0.5, 0.95, 5):
        h = find_hopf(c0, E)
        grid.append(dict(E=float(E), hopf_exists=bool(h),
                         l1=h["l1"] if h else None, supercritical=h["supercritical"] if h else None))
    n_ok = sum(g["hopf_exists"] for g in grid)
    n_super = sum(bool(g["supercritical"]) for g in grid)
    res["robustness"] = {"grid_n": len(grid), "hopf_exists": n_ok, "supercritical": n_super, "grid": grid}
    if base and base["l1"] != 0:
        res["verdict"] = "SUPERCRITICAL" if base["supercritical"] else "SUBCRITICAL"
    else:
        res["verdict"] = "INCONCLUSIVE"
    json.dump(res, open(os.path.join(OUT, "bifurcation_derivation.json"), "w"), indent=2)
    print("== derived Hopf on the frozen fast subsystem (tone control) ==")
    if base:
        print("  Hopf: tau_c=%.4f  G*=%.4f  omega0=%.4f (norm period %.2f)  l1=%.4g -> %s" %
              (base["tau_c"], base["G_c"], base["omega0"], base["period_norm"], base["l1"],
               "SUPERCRITICAL (l1<0)" if base["supercritical"] else "SUBCRITICAL (l1>0)"))
        o = res["onset"]
        print("  onset amplitude^2 vs (tau-tau_c): slope=%.3f R2=%.3f intercept=%.4f -> sqrt-law %s" %
              (o["sqrt_law_slope"], o["sqrt_law_r2"], o["sqrt_law_intercept"],
               "CONSISTENT (supercritical)" if o["consistent_supercritical"] else "not clean"))
    r = res["robustness"]
    print("  robustness across energy settings: Hopf exists %d/%d; supercritical %d/%d; " %
          (r["hopf_exists"], r["grid_n"], r["supercritical"], r["hopf_exists"]))
    print("  VERDICT:", res["verdict"])
    return res


if __name__ == "__main__":
    main()
