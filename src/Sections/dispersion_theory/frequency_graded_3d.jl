function make_frequency_graded_surface(title, table)
    k = table.k_over_sqrt_n0
    x = table.x_over_L
    surface_plot = Plots.surface(k, x, table.lower; color=:blues, colorbar=false,
        title=title, xlabel="k/√n₀", ylabel="x/L", zlabel="ω (rad s⁻¹)",
        xlims=(-pi, pi), ylims=(minimum(x), maximum(x)), zlims=(0, 15),
        camera=(55, 25), size=(1100, 850), dpi=300)
    Plots.surface!(surface_plot, k, x, table.upper; color=:reds, colorbar=false)
    Plots.surface!(surface_plot, k, x, (table.lower .+ table.upper) ./ 2;
        color=:grays, alpha=0.18, colorbar=false)
    return surface_plot
end

function run_frequency_graded_3d(; output_directory=joinpath(
        @__DIR__, "..", "..", "..", "data", "generated", "dispersion_theory",
        "frequency_graded"), n_points=181, n_positions=31)
    mkpath(output_directory)
    zero_tension_plate, resonator = liu_2025_parameters()
    reference_spacing = inv(sqrt(resonator.n0))
    pretension = zero_tension_plate.fluid_density * zero_tension_plate.gravity * reference_spacing
    pretensioned_plate, _ = liu_2025_parameters(tension=pretension)
    grading = FrequencyGradingParameters(resonator.natural_frequency, 1.0, 1.0)
    linear_positions = frequency_graded_positions(grading; n=n_positions,
        profile=ωᵣ_linear)
    exponential_positions = frequency_graded_positions(grading; n=n_positions,
        profile=ωᵣ_exponential)
    normalized_wave_numbers = collect(range(-pi, pi; length=n_points))
    tables = (
        linear_T0=frequency_graded_table(zero_tension_plate, resonator, grading,
            normalized_wave_numbers, linear_positions; profile=ωᵣ_linear),
        linear_Tpositive=frequency_graded_table(pretensioned_plate, resonator, grading,
            normalized_wave_numbers, linear_positions; profile=ωᵣ_linear),
        exponential_T0=frequency_graded_table(zero_tension_plate, resonator, grading,
            normalized_wave_numbers, exponential_positions; profile=ωᵣ_exponential),
        exponential_Tpositive=frequency_graded_table(pretensioned_plate, resonator, grading,
            normalized_wave_numbers, exponential_positions; profile=ωᵣ_exponential),
    )
    titles = (
        linear_T0="Linear grading, T = 0",
        linear_Tpositive="Linear grading, T > 0",
        exponential_T0="Exponential grading, T = 0",
        exponential_Tpositive="Exponential grading, T > 0",
    )
    for name in keys(tables)
        figure = make_frequency_graded_surface(titles[name], tables[name])
        Plots.savefig(figure, joinpath(output_directory,
            "frequency_graded_$(name)_3d.png"))
    end
    return output_directory
end