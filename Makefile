PY ?= python3

.PHONY: verify
verify:
	$(PY) -m pytest -q
	$(PY) -m analysis.nc_bifurcation_derivation
	$(PY) -m analysis.nc_kramers_mfpt
	$(PY) verify_examples.py
