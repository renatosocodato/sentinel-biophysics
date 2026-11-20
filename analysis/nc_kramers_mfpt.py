"""Estimate escape times from the quiet fixed point of the simulated fast subsystem.

Escape is the first crossing of a fixed threshold. Simulations use fixed
seeds and report censoring and bootstrap intervals for observed escape
times. A nominal 240-second cycle converts model time to seconds.
"""
from __future__ import annotations
import json, os
import numpy as np
from scipy.optimize import brentq

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "nc")
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "model"))
from params import P
A, B_, EPS = P["a_fhn"], P["b_fhn"], P["eps_fhn"]
KAPPA = (P["E0"] / (P["K_E"] + P["E0"])) / P["tau_fast"]     # 0.16
TAU = 0.51875                 # fixed tone used by this example
TU_TO_S = 240.0 / 67.34389396762687   # nominal time conversion
G_TH = 0.3                    # escape threshold: focus G*~-0.78, LC peak ~+1 -> 0.3 is a genuine excursion


def focus_fp(tau):
    g = lambda G: G - G**3 / 3 - (G + A) / B_ + tau
    return brentq(g, -1.2, -0.4)                            # the quiescent (left) focus root


def mfpt_focus_escape(sigma, n_ens, seed, Tmax_tu=12000.0, dt=0.05):
    rng = np.random.default_rng(seed)
    Gstar = focus_fp(TAU); wstar = (Gstar + A) / B_
    G = np.full(n_ens, Gstar); w = np.full(n_ens, wstar)
    esc = np.zeros(n_ens, bool); t_esc = np.full(n_ens, np.nan)
    sd = sigma * np.sqrt(dt); nsteps = int(Tmax_tu / dt)
    for i in range(nsteps):
        dG = KAPPA * (G - G**3 / 3 - w + TAU)
        dw = EPS * (G + A - B_ * w)
        G = G + dG * dt + sd * rng.standard_normal(n_ens)
        w = w + dw * dt
        newly = (~esc) & (G > G_TH)
        if newly.any():
            t_esc[newly] = i * dt; esc[newly] = True
            if esc.all():
                break
    return t_esc, esc


def run_sigma(sigma, n_ens=300, seeds=(11, 23, 37)):
    times, ncens = [], 0
    for s in seeds:
        t_esc, esc = mfpt_focus_escape(sigma, n_ens, s)
        times.extend(t_esc[esc].tolist()); ncens += int((~esc).sum())
    times = np.array(times)
    if len(times) < 5:
        return dict(sigma=sigma, mfpt_tu=None, mfpt_s=None, censored_frac=ncens / (n_ens * len(seeds)),
                    note="mostly censored (escape slower than Tmax) -> long-memory / stable focus")
    # bootstrap CI over the pooled escape times
    rng = np.random.default_rng(101)
    boot = [np.mean(rng.choice(times, len(times), replace=True)) for _ in range(400)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return dict(sigma=sigma, n=len(times), censored_frac=ncens / (n_ens * len(seeds)),
                mfpt_tu=float(np.mean(times)), mfpt_s=float(np.mean(times) * TU_TO_S),
                ci95_s=[float(lo * TU_TO_S), float(hi * TU_TO_S)],
                median_s=float(np.median(times) * TU_TO_S))


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {"provenance": "seeded simulation of first-passage times in the frozen model",
           "tau": TAU, "focus_G": focus_fp(TAU), "escape_threshold_G": G_TH, "tu_to_s": TU_TO_S,
           "nominal_cycle_s": 240.0, "by_sigma": {}}
    print("== Simulated first-passage times with bootstrap intervals ==")
    for sigma in (0.02, 0.035, 0.05, 0.08):
        r = run_sigma(sigma)
        res["by_sigma"][f"sigma_{sigma}"] = r
        if r.get("mfpt_s"):
            print("  sigma=%.3f: mean = %.0f s (95%% CI %.0f–%.0f), median %.0f s, censored %.0f%%"
                  % (sigma, r["mfpt_s"], r["ci95_s"][0], r["ci95_s"][1], r["median_s"], 100 * r["censored_frac"]))
        else:
            print("  sigma=%.3f: %s (censored %.0f%%)" % (sigma, r["note"], 100 * r["censored_frac"]))
    with open(os.path.join(OUT, "kramers_mfpt.json"), "w") as output:
        json.dump(res, output, indent=2)
    return res


if __name__ == "__main__":
    main()
