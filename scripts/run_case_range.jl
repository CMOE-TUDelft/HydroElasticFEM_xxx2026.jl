# Compute the simulation records of a range of cases of one test, without
# post-processing, so that heavy tests can be split over resource-limited
# processes. Afterwards `run_tests(id)` only loads the records.
#
#   julia --project -O1 scripts/run_case_range.jl L1-3          # list the cases
#   julia --project -O1 scripts/run_case_range.jl L1-3 1:12     # compute cases 1–12

using HydroElasticFEM_xxx2026

id = ARGS[1]
range = length(ARGS) >= 2 ? eval(Meta.parse(ARGS[2])) : nothing
run_case_range(id, range)
