"""Reference parameters for the six-variable model.

They define numerical examples and are not parameter estimates from
measurements. The unit-interval and dimensional variables are described in
microglia_state_model.py.
"""

P = {
    # ---- Fast contractile oscillator (G, w): FitzHugh-Nagumo scaffold ----
    "tau_fast": 5.0,        # s, fast-variable timescale          [I; from review Box7 tau_fast 2-5 s]
    "eps_fhn": 0.08,        # timescale separation                [I; Melo 2022 EJCB scaffold eps~0.08 (assembled est.)]
    "a_fhn": 0.7,           # FHN offset                          [U; standard FHN prior]
    "b_fhn": 0.8,           # FHN recovery coupling               [U; standard FHN prior]
    "k_c": 0.5,             # GEF-H1 release gain from (1-c)      [U; to be constrained]

    # ---- Tubulin-PTM hysteresis (c): asymmetric write/erase ----
    # tau_c ~ HOURS envelope (Szyk 2014 writer kcat ~1.6/h; Skoge&Ziegler 2017
    # eraser ~1500x slower on polymer): per-second rates must be ~1e-4, not ~1e-2.
    # Values below balance write~erase at c~0.8 in the dark (slow memory).
    "k_write": 3.0e-4,      # /s, alphaTAT1 write rate            [I; Szyk 2014 Cell envelope; abs. microglial rate U]
    "K_write": 0.3,         # write half-saturation               [U]
    "k_erase_base": 1.5e-4, # /s, basal HDAC6 erase               [I; Skoge&Ziegler 2017 SciRep envelope; abs. U]
    "k_erase_gain": 3.0,    # contractile-driven erase boost      [I; depolymerization-coupled erasure (Skoge&Ziegler 2017)]

    # ---- Brake capacity (B): review Box 8 ----
    "k_restore": 0.005,     # /s, brake recovery                  [U]
    "k_desens": 0.02,       # /s, contractile desensitization     [U]
    "K_B": 0.6, "n_B": 3.0, # Hill threshold/cooperativity        [I; review Box8 n=3]

    # ---- Energy charge (E): review Box 3 / 7.13 ----
    # SIGN robust (Socodato 2023 + Cai/Guo 2025 + Sabogal-Guaqueta 2023: activation/RhoA-low
    # -> ATP down, glycolysis up); MAGNITUDE absent (no microglial per-cell flux/EC). E in [0,1]
    # is the Atkinson 1968 energy-charge normalization (0.7-0.95 general prior, not microglial).
    "E_supply": 0.9,        # normalized supply (0..1)            [U; sign-sourced, magnitude unmeasured]
    "e_demand_G": 0.3,      # contractile energy demand           [U]
    "e_demand_c": 0.1,      # PTM-maintenance demand              [U]
    "tau_E": 30.0,          # s, energy relaxation                [U]
    "K_E": 0.2,             # energy gate half-saturation         [I; AMPK feedback form, Muraleedharan&Dasgupta 2022 FEBS J]

    # ---- Process length (L): buckling-limited (review Box4 L_max) ----
    "L_max": 23.0,          # um, max stable length               [L/I; review Box4 ~23 um <-> Biermeier 25um thresh [NEEDS-VERIFY]]
    "k_L": 1.0,             # contractility->shortening gain       [U]
    "tau_L": 60.0,          # s, length relaxation                [I; review tau_slow 60-120 s]

    # ---- Initial conditions (surveillance-like start) ----
    "G0": -1.0, "w0": -0.5, "c0": 0.8, "B0": 0.9, "E0": 0.8, "L0": 20.0,
}
