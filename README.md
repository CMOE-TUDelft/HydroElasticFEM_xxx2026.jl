# HydroElasticFEM_xxx2026.jl

Numerical verification suite for the locally resonant hydroelastic (LRH)
formulation of J2025-LRMM, built on
[HydroElasticFEM.jl](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl).
The tests and their parameters are specified in
[`docs/verification_test_matrix.md`](docs/verification_test_matrix.md).

## Structure

Each test `Lx-y` of the matrix has

- `src/Lx_y_<name>.jl` — computational body: a parameter struct and a
  `run_Lx_y(params)` function returning raw results as a `Dict`;
- `scripts/Lx-<level>/Lx-y-<name>.jl` — case definitions, DrWatson
  bookkeeping, acceptance checks, and the paper figures/tables.

Shared code lives in `src/common/` (workflow utilities, parameter sets,
analytical references, FE builders, observables).

| Location | Content | Committed |
| --- | --- | --- |
| `data/sims/<test-id>/` | raw `produce_or_load` records with git commit/patch and HydroElasticFEM version | no |
| `data/exp_pro/<test-id>/` | processed CSV, including `checks.csv` (pass/fail per acceptance criterion) | yes |
| `plots/<level>/<test-id>/` | paper figures (`.pdf`, `.png`) and tables (`.tex`) | yes |

## Usage

```julia
using HydroElasticFEM_xxx2026
list_tests()                      # registered tests
run_tests("L1-5")                 # one test
run_tests("L1")                   # one level
run_tests("all"; test_suit=true)  # warm-up cases only (smoke test)
```

Run heavy cases with `julia --project -O1 -t <n>`. Tests that register their
case list can be split over several (resource-limited) processes; afterwards
`run_tests(id)` only loads the records and post-processes them:

```bash
julia --project -O1 scripts/run_case_range.jl L1-3          # list the cases
julia --project -O1 scripts/run_case_range.jl L1-3 1:6      # compute cases 1–6
```
 HydroElasticFEM.jl is
tracked from its `main` branch (`Pkg.add(url=..., rev="main")`); every
simulation record stores the Manifest pin (`hefem_rev`, `hefem_tree`).
To work against a local checkout instead, `Pkg.develop(path="../HydroElasticFEM.jl")`.
