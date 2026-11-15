# Scientific scope and limitations

This project demonstrates numerical algorithms and research-software testing with a simplified plasma-inspired model. Important boundaries:

- H₂ level energies can be physically motivated, but the bundled collision tables and most trajectories are explicitly **synthetic**.
- The external Cantera adapter runs a **non-plasma** zero-dimensional NH₃/H₂ reactor if installed successfully. It does not establish coupled plasma-assisted combustion.
- The external two-term Boltzmann EEDF example concerns **N₂/O₂**, not the NH₃/H₂ reactor composition. The EEDF and chemistry adapters are separate.
- Energy-aware graph-channel screening is a **toy proof of concept**, not a validated implementation of the group's plasma-specific P-DRGEP method.
- Two-dimensional transport examples and a separate one-dimensional reactive scalar/thermal proxy are **not** validated reacting CFD/DNS/LES.
- No independent experimental dataset, plasma-assisted ignition data or flame benchmark is supplied. The CSV schema alone is not a validation dataset.
- A mechanism satisfying internal Cantera numerical tolerances would still need independent experimental validation.
- Multiprocessing tests are local workload comparisons; they do not demonstrate HPC strong/weak scaling.

Do not interpret any historical test count or plot as evidence for omitted physics. The work is AI-assisted and the project owner should review every result and reproduce all claimed runs before using it as personal research evidence.
