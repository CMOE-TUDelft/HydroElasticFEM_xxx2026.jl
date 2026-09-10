"Generate every reproduction currently implemented in this repository."
function generate_all_figures(; output_root=joinpath(@__DIR__, "..", "..", "data", "generated"))
    return (non_dissipative_homogeneous_isotropic_lrh_plate=
        run_non_dissipative_homogeneous_isotropic_lrh_plate(output_directory=joinpath(
            output_root, "dispersion_theory", "homogeneous_isotropic")),)
end