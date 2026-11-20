"""Structural invariants of the frozen SENTINEL model.

These are self-contained unit checks on model/ (no external data, no results/ artifacts): they guard
the frozen definition and its declared parameters, so an accidental edit to the model or params is
caught immediately. The numerical examples have their own reference checks in verify_examples.py.
"""
import numpy as np

from model.microglia_state_model import rhs, integrate, observable_complexity, opto_drive
from model.params import P


X0 = [P["G0"], P["w0"], P["c0"], P["B0"], P["E0"], P["L0"]]


def test_state_is_six_dimensional():
    d = rhs(0.0, X0)
    assert np.shape(d) == (6,), "the frozen model is a six-variable system (Gamma, w, c, B, E, L)"


def test_rhs_finite_on_baseline_state():
    assert np.all(np.isfinite(rhs(0.0, X0)))


def test_frozen_baseline_parameters():
    # exactly the declared model values; guards against silent drift
    expected = {
        "tau_fast": 5.0, "eps_fhn": 0.08, "a_fhn": 0.7, "b_fhn": 0.8, "k_c": 0.5,
        "K_write": 0.3, "k_write": 3.0e-4, "k_erase_base": 1.5e-4,
        "K_B": 0.6, "n_B": 3.0, "K_E": 0.2, "E_supply": 0.9, "tau_E": 30.0,
        "L_max": 23.0, "k_L": 1.0, "tau_L": 60.0,
    }
    for k, v in expected.items():
        assert P[k] == v, f"frozen param {k} drifted: {P[k]!r} != {v!r}"


def test_integration_succeeds_and_stays_finite():
    sol = integrate(t_span=(0.0, 240.0))
    assert sol.success
    assert np.all(np.isfinite(sol.y))


def test_observable_in_unit_interval():
    sol = integrate(t_span=(0.0, 240.0))
    for xt in sol.y.T[::10]:
        y = float(observable_complexity(xt))
        assert 0.0 <= y <= 1.0, f"morphological observable out of [0,1]: {y}"


def test_energy_charge_and_length_stay_physical():
    sol = integrate(t_span=(0.0, 600.0))
    E = sol.y[4]
    L = sol.y[5]
    assert np.all(E >= -1e-6) and np.all(E <= 1.0 + 1e-6), "energy charge E must stay in [0,1]"
    assert np.all(L >= -1e-6), "process length L must stay non-negative"


def test_opto_drive_is_a_gated_square_wave():
    # off by default, on within the duty fraction of a period when amplitude>0
    assert opto_drive(0.0) == 0.0
    on = opto_drive(0.0, amplitude=0.4, period=240.0, duty=0.5)
    off = opto_drive(180.0, amplitude=0.4, period=240.0, duty=0.5)
    assert on == 0.4 and off == 0.0
