#!/usr/bin/env bash
# Real chemistry requires an internet-connected Python environment with Cantera.
set -euo pipefail
cd "$(dirname "$0")"
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test,cantera]'
python -c 'import cantera; print("Cantera",cantera.__version__)'
python -m pytest -q tests/test_v11_extensions.py tests/test_v12_regressions.py tests/test_v11_optional_integrations.py
python -m plasmarom ammonia-h2 --out results_v11/ammonia_h2
python -m plasmarom boltzmann --out results_v11/cantera_boltzmann
python -m plasmarom transport-1d --out results_v11/transport
printf '\nComputations completed. Experimental PAC validation and DNS/LES are not claimed.\n'
