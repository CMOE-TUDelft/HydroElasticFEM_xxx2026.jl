module HydroElasticFEM_xxx2026

using Plots
using LaTeXStrings

include("Models/Configurations.jl")
include("Physics/Dispersion.jl")
include("PostProcessing/Dispersion.jl")
include("Utilities/IO.jl")
include("Figures/Figure03_WetModes.jl")
include("Cases/Case01_WetModes.jl")
include("Cases/RunAll.jl")

export PlateParameters, ResonatorParameters, effective_mass,
    bare_dispersion, metaplate_dispersion, liu_2025_parameters,
    run_case01, run_liu_2025_figure3, dispersion_table, bandgap_limits,
    write_dispersion_csv, generate_all_figures

end