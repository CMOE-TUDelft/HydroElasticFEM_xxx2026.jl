"Generate every reproduction currently implemented in this repository."
function generate_all_figures(; output_root=joinpath(@__DIR__, "..", "..", "data", "generated"))
    return (non_dissipative_homogeneous_isotropic_lrh_plate=
        run_non_dissipative_homogeneous_isotropic_lrh_plate(output_directory=joinpath(
            output_root, "dispersion_theory", "homogeneous_isotropic")),
        frequency_graded=run_frequency_graded_dispersion(output_directory=joinpath(
            output_root, "dispersion_theory", "frequency_graded")),
        frequency_graded_3d=run_frequency_graded_3d(output_directory=joinpath(
            output_root, "dispersion_theory", "frequency_graded")),
        density_graded=run_density_graded_dispersion(output_directory=joinpath(
            output_root, "dispersion_theory", "density_graded")),
        density_graded_3d=run_density_graded_3d(output_directory=joinpath(
            output_root, "dispersion_theory", "density_graded")),
        anisotropic_density=run_anisotropic_density_dispersion(output_directory=joinpath(
            output_root, "dispersion_theory", "anisotropic_density")))
end