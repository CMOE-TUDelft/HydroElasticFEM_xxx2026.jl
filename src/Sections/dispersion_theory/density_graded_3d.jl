function make_density_graded_surface(title, table)
    return make_frequency_graded_surface(title, table)
end

function run_density_graded_3d(; output_directory=joinpath(
        @__DIR__, "..", "..", "..", "data", "generated", "dispersion_theory",
        "density_graded"), n_points=181, n_positions=31)
    mkpath(output_directory)
    zero_tension_plate, resonator = liu_2025_parameters()
    reference_spacing = inv(sqrt(resonator.n0))
    pretension = zero_tension_plate.fluid_density * zero_tension_plate.gravity * reference_spacing
    pretensioned_plate, _ = liu_2025_parameters(tension=pretension)
    grading = DensityGradingParameters(resonator.n0, 1.0)
    normalized_positions = density_graded_positions(n=n_positions)
    normalized_wave_numbers = collect(range(-pi, pi; length=n_points))
    tables = (
        T0=density_graded_table(zero_tension_plate, resonator, grading,
            normalized_wave_numbers, normalized_positions),
        Tpositive=density_graded_table(pretensioned_plate, resonator, grading,
            normalized_wave_numbers, normalized_positions),
    )
    titles = (T0="Density grading, T = 0", Tpositive="Density grading, T > 0")
    for name in keys(tables)
        figure = make_density_graded_surface(titles[name], tables[name])
        PlotlyJS.savefig(figure, joinpath(output_directory,
            "density_graded_$(name)_3d.html"))
    end
    return output_directory
end