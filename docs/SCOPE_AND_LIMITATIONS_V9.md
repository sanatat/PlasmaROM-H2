# Scope and evidence boundary (v9)

This repository is an **educational, AI-assisted numerical chemistry research prototype**. Local verification covers ODE/PDE numerical behavior, synthetic collision-rate integrations, scenario separation, and mock-tested acceptance invariants. It is *not* an experimentally validated plasma-assisted-combustion solver.

The optional real-chemistry pipeline uses Cantera, but has not been run locally: there is no accepted H2/O2/NH3 reduced mechanism in this release. Publication-worthy validation requires independently documented electron-impact rates, EEDF/Boltzmann solution, state-to-state kinetic mechanism, energy-partition-preserving P-DRGEP, ignition/flame validation and convergence/error analysis. DNS/LES and MPI/GPU solver scaling also remain outside scope.

For a graduate application, report only modules personally run, inspected and understood. Do not claim experimental plasma experience, P-DRGEP validation, true turbulent DNS/LES, or HPC acceleration from these synthetic demos.
