function run_frequency_graded_dispersion(; output_directory=joinpath(
        @__DIR__, "..", "..", "..", "data", "generated", "dispersion_theory",
        "frequency_graded"), n_points=701, n_positions=7)
    gr()
    mkpath(output_directory)
    normalized_wave_numbers = collect(range(-pi, pi; length=n_points))
    zero_tension_plate, resonator = liu_2025_parameters()
    reference_spacing = inv(sqrt(resonator.n0))
    pretension = zero_tension_plate.fluid_density * zero_tension_plate.gravity * reference_spacing
    pretensioned_plate, _ = liu_2025_parameters(tension=pretension)
    grading = FrequencyGradingParameters(resonator.natural_frequency, 1.0, 1.0)
    normalized_positions = frequency_graded_positions(grading; n=n_positions,
        profile=ωᵣ_linear)
    linear = (plate -> frequency_graded_table(plate, resonator, grading,
        normalized_wave_numbers, normalized_positions; profile=ωᵣ_linear))
    exponential_positions = frequency_graded_positions(grading; n=n_positions,
        profile=ωᵣ_exponential)
    exponential = (plate -> frequency_graded_table(plate, resonator, grading,
        normalized_wave_numbers, exponential_positions; profile=ωᵣ_exponential))
    tables = (linear_zero=linear(zero_tension_plate),
        linear_positive=linear(pretensioned_plate),
        exponential_zero=exponential(zero_tension_plate),
        exponential_positive=exponential(pretensioned_plate))
    plots = (
        linear_T0=make_frequency_graded_panel(L"\mathrm{Linear\ grading},\ T=0",
            tables.linear_zero),
        linear_Tpositive=make_frequency_graded_panel(
            L"\mathrm{Linear\ grading},\ T>0", tables.linear_positive),
        exponential_T0=make_frequency_graded_panel(
            L"\mathrm{Exponential\ grading},\ T=0", tables.exponential_zero),
        exponential_Tpositive=make_frequency_graded_panel(
            L"\mathrm{Exponential\ grading},\ T>0", tables.exponential_positive),
    )
    for (name, figure) in pairs(plots)
        savefig(figure, joinpath(output_directory, "frequency_graded_$(name).png"))
    end
    return output_directory
end