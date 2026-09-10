module HydroElasticFEM_xxx2026

using Plots
using LaTeXStrings

include("Models/Configurations.jl")
include("Physics/Dispersion.jl")
include("PostProcessing/Dispersion.jl")
include("Utilities/IO.jl")
include("Sections/dispersion_theory/homogeneous_isotropic.jl")
include("Sections/dispersion_theory/homogeneous_isotropic_plot.jl")
include("Cases/RunAll.jl")

export PlateParameters, ResonatorParameters, effective_mass,
    bare_dispersion, metaplate_dispersion, liu_2025_parameters,
    run_non_dissipative_homogeneous_isotropic_lrh_plate, dispersion_table, bandgap_limits,
    write_dispersion_csv, generate_all_figures

end