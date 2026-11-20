"""Reference parameters for the six-variable model.

They define numerical examples and are not parameter estimates from
measurements. The unit-interval and dimensional variables are described in
microglia_state_model.py.
"""

P = {
    # ---- Fast contractile oscillator (G, w): FitzHugh-Nagumo scaffold ----
    "tau_fast": 5.0,        # s, fast-variable timescale
    "eps_fhn": 0.08,        # timescale separation
    "a_fhn": 0.7,           # FHN offset
    "b_fhn": 0.8,           # FHN recovery coupling
    "k_c": 0.5,             # GEF-H1 release gain from (1-c)

    # ---- Tubulin-PTM hysteresis (c): asymmetric write/erase ----
    "k_write": 3.0e-4,      # /s, alphaTAT1 write rate
    "K_write": 0.3,         # write half-saturation
    "k_erase_base": 1.5e-4, # /s, basal HDAC6 erase
    "k_erase_gain": 3.0,    # contractile-driven erase boost

    # ---- Brake capacity (B)----
    "k_restore": 0.005,     # /s, brake recovery
    "k_desens": 0.02,       # /s, contractile desensitization
    "K_B": 0.6, "n_B": 3.0, # Hill threshold/cooperativity

    # ---- Energy charge (E)----
    "E_supply": 0.9,        # normalized supply (0..1)
    "e_demand_G": 0.3,      # contractile energy demand
    "e_demand_c": 0.1,      # PTM-maintenance demand
    "tau_E": 30.0,          # s, energy relaxation
    "K_E": 0.2,             # energy gate half-saturation

    # ---- Process length (L): buckling-limited ----
    "L_max": 23.0,          # um, max stable length
    "k_L": 1.0,             # contractility->shortening gain
    "tau_L": 60.0,          # s, length relaxation

    # ---- Initial conditions (surveillance-like start) ----
    "G0": -1.0, "w0": -0.5, "c0": 0.8, "B0": 0.9, "E0": 0.8, "L0": 20.0,
}
