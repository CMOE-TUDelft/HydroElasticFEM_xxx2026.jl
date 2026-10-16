"""
Numerical verification of the locally resonant hydroelastic (LRH) formulation
of the paper J2025-LRMM, built on HydroElasticFEM.jl.

Layout (see README.md):

- `src/references.jl` — the only shared code: analytical references of the
  paper (Sections 3–4), no FE;
- `src/Lx_y_<name>.jl` — one self-contained module per test: its parameter
  struct (fixed and variable parameters), the FE problem written with the
  HydroElasticFEM.jl API, and the observables of one case;
- `scripts/<level>/Lx-y_<name>.jl` — one script per test: the cases of the
  study, DrWatson records, post-processing, checks, figures and tables.

Tests are specified in `docs/verification_test_matrix.md`.
"""
module HydroElasticFEM_xxx2026

include("references.jl")
include("L1_1_bare_host_dispersion.jl")
include("L1_2_flux_normalisation.jl")       # uses the L1-1 problem
include("L1_3_single_line_resonator.jl")    # L1-1 problem with a resonator

export References, L1_1, L1_2, L1_3

end
