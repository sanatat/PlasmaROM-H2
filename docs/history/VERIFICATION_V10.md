# Version 10 verification — 2026-11-15

**Status:** unit and adversarial behavior checked locally; optional physical chemistry remains unverified without Cantera.

## Local execution

| Suite | Pass | Skip | Failure | Notes |
|---|---:|---:|---:|---|
| Non-scaling pytest | 140 | 5 | 0 | Cantera unavailable; see `results_v10/PYTEST_CORE.log` |
| Scaling module, separate | 2 | 0 | 0 | See `results_v10/PYTEST_SCALING.log` |
| Targeted v10 adversarial tests | 11 | 0 | 0 | Subset of non-scaling suite |

**Not claimed:** an uninterrupted single full-suite run; the combined run hit a multiprocessing-timeout boundary. Two separate successful runs are reported, not added into a fictional combined log.

## Verified v10 behavior

- Collider efficiency and explicit collision partners change gross flux connectivity; generic M defaults are handled conservatively.
- A tiny OH reference `1e-12` and larger candidate `5e-9` is now flagged by the trace relative metric, despite the previous normalised error staying below 0.12.
- Training replay altered in a mocked Cantera layer is rejected, even when original JSON aggregates remain internally consistent.
- Missing training candidate records and unpinned protocol markers fail the v10 integrity gate.
- Pinned evaluation plans accept matching run conditions and reject altered tolerance / observation horizon.

## Scientific limitations

Mocked or synthetic tests verify **logic**, not the precision of real Cantera integration. The chemistry tests remain skipped because Cantera is absent; pip cannot access its package index from this environment. No chemically validated reduction YAML is included. No external electron-cross-section datasets are included. A failing scientific gate is an honest negative result.
