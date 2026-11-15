# PlasmaROM-H2 v11 — evidence status and unmet research goals

This document separates **written**, **tested locally**, **tested only with mocks**, and **not scientifically demonstrated**. The project is an AI-assisted prototype, not an experimental or validated PAC simulation.

| Task | Current state | What evidence exists / is missing |
|---|---|---|
| Cantera real chemical kinetics | IMPLEMENTED, NOT RUN LOCALLY | New NH3/H2 reactor uses documented Alzueta-2023 example data (phase `baseline`); Cantera unavailable offline. GitHub Actions requires Cantera to import and execute. |
| Two-term Boltzmann EEDF | IMPLEMENTED AS CANTERA 3.2 ADAPTER, NOT RUN LOCALLY | Official Phelps N2/O2 example + 50, 100, 200 Td; **not H2/NH3**; checks EEDF normalization. No discharge self-consistency. |
| Real electron collision dataset | AVAILABLE VIA EXTERNAL CANTERA EXAMPLE; NOT BUNDLED | Phelps compilation and Itikawa 2009 examples are traceable, subject to source-specific attribution. Synthetic H2 data retained in legacy examples. |
| Plasma P-DRGEP | NOT COMPLETED | New energy-aware graph screen verified on synthetic data; **not** the published algorithm and no real plasma state-network validation. |
| NH3/H2 reactor numerical validation | IMPLEMENTED, NOT RUN LOCALLY | Alzueta 2023 reactor code, ignition censoring, species outputs, provenance. No experimental result or accepted reduction. |
| Reactive flow | 1D SCALAR TRANSPORT EXECUTED | Conservative species/energy transport with analytical homogeneous-reaction check. **Not** variable-density reacting Navier–Stokes or combustion CFD. |
| DNS/LES | NOT COMPLETED | Neither 3D turbulence nor validated LES sub-grid closure. Existing 2D flow remains educational. |
| Experimental PAC validation | NOT COMPLETED | Added a strict CSV measurement/model comparison tool. No measured data, uncertainties, or calibration from an independent PAC experiment included. |
| HPC | NOT COMPLETED | Existing repeated local multi-processing benchmarks; no solver scaling on HPC cluster. |

## External sources (data not redistributed)

- Cantera Phelps two-term tutorial: https://cantera.org/stable/examples/python/thermo/plasma-eedf.html
- Cantera ammonia model: https://github.com/Cantera/cantera-example-data/blob/main/ammonia-CO-H2-Alzueta-2023.yaml
- Cantera Itikawa oxygen example: https://cantera.org/stable/examples/input/oxygen-plasma-itikawa.html
- BOLOS solver project: https://github.com/aluque/bolos
- LXCat data provenance/rights: https://master.lxcat.net/instructions/how_reference.php

## Academic honesty for applications

Appropriate: **Numerical verification, conservative 1D reactive transport, graph/energy screening prototyping, and implementation of Cantera/Boltzmann integration workflows.**

Not appropriate: claiming **completed experimentally validated NH3/H2 plasma combustion, published P-DRGEP, 3D DNS/LES, successful mechanism reduction or observed speedup**. A smoke test is not experimental validation.
