"""Check the selected reference values for the synthetic model examples."""

import json
import math
from pathlib import Path


def main() -> None:
    output = Path(__file__).resolve().parent / "results" / "nc"
    hopf = json.loads((output / "bifurcation_derivation.json").read_text())["hopf"]
    escape = json.loads((output / "kramers_mfpt.json").read_text())["by_sigma"]
    cases = (
        ("Hopf tone", hopf["tau_c"], 0.5265, 0.001),
        ("First Lyapunov coefficient", hopf["l1"], 0.2355, 0.01),
        ("Mean escape time at noise 0.02 (s)", escape["sigma_0.02"]["mfpt_s"], 2066.55, 5.0),
        ("Mean escape time at noise 0.05 (s)", escape["sigma_0.05"]["mfpt_s"], 510.51, 5.0),
    )
    for label, actual, reference, tolerance in cases:
        if not math.isclose(actual, reference, rel_tol=0, abs_tol=tolerance):
            raise SystemExit(f"{label}: {actual:g}, expected {reference:g} ± {tolerance:g}")
        print(f"{label}: {actual:g} (reference check passed)")


if __name__ == "__main__":
    main()
