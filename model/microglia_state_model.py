"""
microglia_state_model.py — coupled slow-fast model of microglial state (companion paper).

STATUS: UNCALIBRATED v0 (Milestone 0). This module SPECIFIES and INTEGRATES the model;
it does NOT validate it. Parameter values are placeholders/priors sourced in
PARAMETER_PROVENANCE.md (many flagged NEEDS-VERIFY / U=unconstrained). The central
hypotheses are to be TESTED in later milestones, not assumed here:
  H1: the coupled system supports surveillance(limit cycle)/reactive(fixed pt)/
      dystrophic(collapse)/quiescent(node) as coexisting attractors  -> NOT YET TESTED
  H2: periodic opto-RhoA forcing entrains the surveillance oscillator -> NOT YET TESTED
A negative result for H1/H2 is an acceptable, reportable outcome.

State vector x = (G, w, c, B, E, L):
  G  = Gamma, contractile balance (fast activator; opto-RhoA acts here)   [dimensionless, FHN-scaled]
  w  = fast recovery/inhibitor (FitzHugh-Nagumo recovery variable)        [dimensionless]
  c  = tubulin-PTM (acetylation) landscape (slow memory; 0..1)            [fraction acetylated]
  B  = brake capacity (slow adaptation; 0..1)                             [fraction responsive]
  E  = energy charge epsilon_E (gate; 0..1)                               [normalized ATP availability]
  L  = process length (mechanical, slow)                                  [um]

Conceptual mapping to the review (frozen at Stage 2): G<->Box4 Gamma & Box7 activator u;
w<->Box7 inhibitor v; c<->Box6 tubulin-code hysteresis; B<->Box8 brake; E<->Box3/7.13 energy;
L<->Box4 L_max / 6.1.2 internalization threshold. opto-RhoA = drive I(t) (Box 10).

Requires: numpy, scipy.
"""
from __future__ import annotations
import numpy as np

try:
    from .params import P  # parameter dict with provenance tags
except Exception:  # allow running as a script without package context
    from params import P


# ----------------------------------------------------------------------------
# Forcing: opto-RhoA actuator I(t). This is the falsifiability hook (Box 10).
# ----------------------------------------------------------------------------
def opto_drive(t, amplitude=0.0, period=None, duty=0.5, t_on=0.0, t_off=np.inf):
    """opto-ArhGEF25 drive into the fast (Gamma) equation.
    amplitude=0 -> dark (no stimulation). period=None -> constant clamp while on.
    For the ENTRAINMENT test (H2, later milestone): sweep `period` across the
    intrinsic surveillance band and measure phase-locking. Here it is only a hook.
    """
    if amplitude == 0.0 or t < t_on or t > t_off:
        return 0.0
    if period is None:
        return amplitude
    phase = ((t - t_on) % period) / period
    return amplitude if phase < duty else 0.0


# ----------------------------------------------------------------------------
# Right-hand side of the ODE system.
# ----------------------------------------------------------------------------
def rhs(t, x, p=P, drive_kwargs=None):
    G, w, c, B, E, L = x
    dk = drive_kwargs or {}
    I = opto_drive(t, **dk)

    # --- Fast layer: FitzHugh-Nagumo-type contractile oscillator (G, w) ---
    # Low acetylation (1-c) releases GEF-H1 -> biases Gamma up (review Box 5/6).
    # Energy gates the achievable contractile output (review Box 3 AMPK cap).
    gef_release = p["k_c"] * (1.0 - c)
    energy_gate = E / (p["K_E"] + E)
    dG = (G - G**3 / 3.0 - w + I + gef_release) * energy_gate / p["tau_fast"]
    dw = p["eps_fhn"] * (G + p["a_fhn"] - p["b_fhn"] * w)

    # --- Intermediate slow memory: tubulin-PTM hysteresis (c) ---
    # Asymmetric write/erase (review Box 6): slow lumen-limited writing,
    # fast contractile-activation-driven erasing.
    write = p["k_write"] * (1.0 - c) / (p["K_write"] + (1.0 - c))
    erase = p["k_erase_base"] * (1.0 + p["k_erase_gain"] * max(G, 0.0)) * c
    dc = write - erase

    # --- Slow adaptation: brake capacity (B), review Box 8 ---
    desens = p["k_desens"] * (max(G, 0.0) ** p["n_B"]) / (p["K_B"] ** p["n_B"] + max(G, 0.0) ** p["n_B"])
    dB = p["k_restore"] * (1.0 - B) - desens * B

    # --- Energy charge (E in [0,1]), review Box 3 / 7.13; TREM2-mTOR supply ---
    # Logistic supply/demand keeps E bounded in [0,1] for E_supply in [0,1]:
    # supply fills the deficit (1-E); cytoskeletal demand drains proportional to E.
    demand = p["e_demand_G"] * max(G, 0.0) + p["e_demand_c"] * (1.0 - c)
    dE = (p["E_supply"] * (1.0 - E) - demand * E) / p["tau_E"]

    # --- Process length (L), buckling-limited (review Box4 L_max / 6.1.2) ---
    # Higher contractility -> shorter processes; capped by L_max.
    L_target = p["L_max"] / (1.0 + p["k_L"] * max(G, 0.0))
    dL = (L_target - L) / p["tau_L"]

    return [dG, dw, dc, dB, dE, dL]


def observable_complexity(x, p=P):
    """Map state -> measurable morphological order parameter (Biermeier 2026):
    complexity ~ 1/solidity. High Gamma (contractile) + short L -> low complexity.
    Returned in [0,1] (1 = high-complexity ramified, 0 = compact)."""
    G, w, c, B, E, L = (x if np.ndim(x) == 1 else x)
    g = np.clip((np.asarray(G) + 2.0) / 4.0, 0, 1)         # normalize FHN range ~[-2,2]
    lnorm = np.clip(np.asarray(L) / p["L_max"], 0, 1)
    return np.clip(0.5 * (1.0 - g) + 0.5 * lnorm, 0, 1)


def integrate(t_span=(0.0, 1800.0), x0=None, p=P, drive_kwargs=None, max_step=1.0):
    """Integrate the system. Returns scipy solution object.
    NOTE: this only demonstrates that the system integrates; it makes NO claim
    about attractors or entrainment (those are later-milestone analyses)."""
    from scipy.integrate import solve_ivp
    if x0 is None:
        x0 = [p["G0"], p["w0"], p["c0"], p["B0"], p["E0"], p["L0"]]
    return solve_ivp(rhs, t_span, x0, args=(p, drive_kwargs),
                     method="LSODA", max_step=max_step, dense_output=True, rtol=1e-6, atol=1e-9)


if __name__ == "__main__":
    # ---- Milestone-0 structural sanity check ONLY (no scientific claim) ----
    sol = integrate()
    ok = bool(sol.success)
    xf = sol.y[:, -1]
    comp = float(observable_complexity(xf))
    print("[microglia_state_model v0] integration success:", ok)
    print("  final state (G,w,c,B,E,L):", np.round(xf, 4).tolist())
    print("  final complexity proxy (0=compact,1=ramified):", round(comp, 4))
    print("  NOTE: uncalibrated v0; this is a structural sanity check, NOT a result.")
    print("  H1 (attractor coexistence) and H2 (entrainment) are UNTESTED here.")
