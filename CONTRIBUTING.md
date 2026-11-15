# Contributing

Contributions are welcome when they strengthen reproducibility, numerical correctness or documented physical scope.

1. Open an issue describing the proposed change, the equations or evidence, and a minimal reproducible example.
2. Run `python -m pytest -q tests/test_v12_regressions.py` and applicable targeted tests. Supply an independent reference for physical claims.
3. New datasets must include source, DOI/URL, units, uncertainties, usage license, processing steps, and whether measurements are independent. Keep synthetic and experimental data clearly separated.
4. For new numerical results, specify the parameter ranges, solver tolerances, grid/time convergence, test/holdout separation and any failures.
5. Never claim that green CI, small residuals or synthetic validation establish experimental plasma-assisted combustion or DNS/LES fidelity.

Optional Cantera, Boltzmann and experimental integrations require machine-readable artifacts and source citations. Do not commit secrets, API tokens or proprietary research data.
