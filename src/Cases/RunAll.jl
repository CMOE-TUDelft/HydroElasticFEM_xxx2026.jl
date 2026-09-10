"Generate every reproduction currently implemented in this repository."
function generate_all_figures(; output_root=joinpath(@__DIR__, "..", "..", "data", "generated"))
    return (case01=run_case01(output_directory=joinpath(output_root,
        "subsubsec_non_dissipative_homogeneous_isotropic_dispersion")),)
end