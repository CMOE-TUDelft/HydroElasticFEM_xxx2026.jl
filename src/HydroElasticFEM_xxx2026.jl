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
include("Cases/RunAll.jl")

export PlateParameters, ResonatorParameters, FrequencyGradingParameters,
    effective_mass,
    bare_dispersion, metaplate_dispersion, liu_2025_parameters,
    run_non_dissipative_homogeneous_isotropic_lrh_plate, dispersion_table, bandgap_limits,
    write_dispersion_csv, ωᵣ_linear, ωᵣ_exponential, graded_dispersion,
    frequency_graded_table, run_frequency_graded_dispersion,
    make_frequency_graded_surface, run_frequency_graded_3d, generate_all_figures

end