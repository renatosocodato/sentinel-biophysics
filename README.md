# SENTINEL

SENTINEL is a computational model of microglial surveillance. It links a fast
contractile oscillator to slower changes in the cytoskeleton, energy and process
length. This release provides the equations, their structural tests and two
self-contained numerical examples.

## Reproduce the examples

With Python 3.14 and Make installed, run these commands from the repository root.

```sh
python3.14 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 make verify PY=.venv/bin/python
```

The command runs the structural tests and regenerates two JSON files in
`results/nc/`. The first example locates a Hopf crossing and checks its type
using the first Lyapunov coefficient. The second estimates mean first-passage
times in seeded, model-generated trajectories at different noise amplitudes,
mapping a linear reference period of about 67.34 model time units to a
nominal 240 seconds. The escape-time fields ending in `_s` use that scale.
The verifier compares four values with the reference values and absolute
tolerances written in `verify_examples.py`.

The examples characterize the specified equations and parameter settings.
Their results are numerical model outputs. Both run without experimental data.

## Cite this software

Socodato R, Relvas JB (2026). *SENTINEL: computational biophysics of microglial
surveillance*, v0.1.0. See `CITATION.cff` for the version-specific citation.

The software is distributed under the MIT license.
