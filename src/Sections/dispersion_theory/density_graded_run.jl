function run_density_graded_dispersion(; output_directory=joinpath(
        @__DIR__, "..", "..", "..", "data", "generated", "dispersion_theory",
        "density_graded"), n_points=701, n_positions=7)
    Plots.gr()
    mkpath(output_directory)
    zero_tension_plate, resonator = liu_2025_parameters()
    reference_spacing = inv(sqrt(resonator.n0))
    pretension = zero_tension_plate.fluid_density * zero_tension_plate.gravity * reference_spacing
    pretensioned_plate, _ = liu_2025_parameters(tension=pretension)
    grading = DensityGradingParameters(resonator.n0, 1.0)
    normalized_wave_numbers = collect(range(-pi, pi; length=n_points))
    normalized_positions = density_graded_positions(n=n_positions)
    tables = (T0=density_graded_table(zero_tension_plate, resonator, grading,
        normalized_wave_numbers, normalized_positions),
        Tpositive=density_graded_table(pretensioned_plate, resonator, grading,
            normalized_wave_numbers, normalized_positions))
    plots = (T0=make_density_graded_panel(L"\mathrm{Density\ grading},\ T=0",
            tables.T0),
        Tpositive=make_density_graded_panel(L"\mathrm{Density\ grading},\ T>0",
            tables.Tpositive))
    for (name, figure) in pairs(plots)
        Plots.savefig(figure, joinpath(output_directory, "density_graded_$(name).png"))
    end
    return output_directory
end