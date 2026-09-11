module HydroElasticFEM_xxx2026

using Plots
using LaTeXStrings
using PlotlyJS

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
    anisotropic_optical_depths, run_anisotropic_density_dispersion,
    generate_all_figures

end