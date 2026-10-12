module HydroElasticFEM_xxx2026

using DrWatson
using CairoMakie
using DataFrames
using CSV
using Printf

# ─── Legacy Section-3 local-theory figures (Plots.jl) ───────────────────────
# Kept in their own module so that Plots.jl and CairoMakie do not clash.
# To be migrated to the S3-x test pattern (verification test matrix, §7).
module LegacyDispersionTheory

using Plots
using LaTeXStrings

include("Models/Configurations.jl")
include("Physics/Dispersion.jl")
include("PostProcessing/Dispersion.jl")
include("Utilities/IO.jl")
include("Sections/dispersion_theory/homogeneous_isotropic.jl")
include("Sections/dispersion_theory/homogeneous_isotropic_plot.jl")
include("Sections/dispersion_theory/frequency_graded.jl")
include("Sections/dispersion_theory/frequency_graded_plot.jl")
include("Sections/dispersion_theory/frequency_graded_run.jl")
include("Sections/dispersion_theory/frequency_graded_3d.jl")
include("Sections/dispersion_theory/density_graded.jl")
include("Sections/dispersion_theory/density_graded_plot.jl")
include("Sections/dispersion_theory/density_graded_run.jl")
include("Sections/dispersion_theory/density_graded_3d.jl")
include("Sections/dispersion_theory/anisotropic_density.jl")
include("Sections/dispersion_theory/anisotropic_density_plot.jl")
include("Sections/dispersion_theory/anisotropic_density_run.jl")
include("Cases/RunAll.jl")

export PlateParameters, ResonatorParameters, FrequencyGradingParameters,
    DensityGradingParameters,
    effective_mass,
    bare_dispersion, metaplate_dispersion, liu_2025_parameters,
    run_non_dissipative_homogeneous_isotropic_lrh_plate, dispersion_table, bandgap_limits,
    write_dispersion_csv, ωᵣ_linear, ωᵣ_exponential, graded_dispersion,
    frequency_graded_table, run_frequency_graded_dispersion,
    make_frequency_graded_surface, run_frequency_graded_3d,
    n_density_graded, Mᵣ, density_graded_dispersion, density_graded_table,
    run_density_graded_dispersion, make_density_graded_surface, run_density_graded_3d,
    anisotropic_density, anisotropic_resonant_wavenumber, anisotropic_density_field,
    anisotropic_density_table, make_anisotropic_density_panel,
    anisotropic_dispersion_field, make_anisotropic_density_surface,
    anisotropic_optical_depths,
    run_anisotropic_density_dispersion,
    generate_all_figures

end # module LegacyDispersionTheory

using .LegacyDispersionTheory

export PlateParameters, ResonatorParameters, FrequencyGradingParameters,
    DensityGradingParameters,
    effective_mass,
    bare_dispersion, metaplate_dispersion, liu_2025_parameters,
    run_non_dissipative_homogeneous_isotropic_lrh_plate, dispersion_table, bandgap_limits,
    write_dispersion_csv, ωᵣ_linear, ωᵣ_exponential, graded_dispersion,
    frequency_graded_table, run_frequency_graded_dispersion,
    make_frequency_graded_surface, run_frequency_graded_3d,
    n_density_graded, Mᵣ, density_graded_dispersion, density_graded_table,
    run_density_graded_dispersion, make_density_graded_surface, run_density_graded_3d,
    anisotropic_density, anisotropic_resonant_wavenumber, anisotropic_density_field,
    anisotropic_density_table, make_anisotropic_density_panel,
    anisotropic_dispersion_field, make_anisotropic_density_surface,
    anisotropic_optical_depths,
    run_anisotropic_density_dispersion,
    generate_all_figures

# ─── Verification suite (docs/verification_test_matrix.md) ──────────────────
# Each test Lx-y has a computational body in `src/Lx_y_<name>.jl` and a
# pre/post-processing script in `scripts/Lx-<level>/Lx-y-<name>.jl` that
# registers its entry point with `register_test`.

include("plot_theme.jl")
include("common/workflow.jl")

export run_tests, list_tests
export output_root, with_output_root, sims_path, exp_pro_path, plots_path
export provenance, produce_or_load_case
export AcceptanceCheck, add_check!, write_checks, save_figure, write_latex_table

end
