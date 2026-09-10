function make_frequency_graded_surface(title, table)
    k = table.k_over_sqrt_n0
    x = table.x_over_L
    lower_surface = PlotlyJS.surface(x=k, y=x, z=table.lower,
        colorscale=[[0.0, "royalblue"], [1.0, "royalblue"]],
        showscale=false, opacity=0.9, name="lower branch",
        hovertemplate="k/√n₀=%{x:.3f}<br>x/L=%{y:.3f}<br>ω₋=%{z:.3f}<extra></extra>")
    upper_surface = PlotlyJS.surface(x=k, y=x, z=table.upper,
        colorscale=[[0.0, "firebrick"], [1.0, "firebrick"]],
        showscale=false, opacity=0.9, name="upper branch",
        hovertemplate="k/√n₀=%{x:.3f}<br>x/L=%{y:.3f}<br>ω₊=%{z:.3f}<extra></extra>")
    bandgap_surface = PlotlyJS.surface(x=k, y=x,
        z=(table.lower .+ table.upper) ./ 2,
        colorscale=[[0.0, "rgba(120,120,120,0.18)"],
            [1.0, "rgba(120,120,120,0.18)" ]], showscale=false, opacity=0.18,
        name="bandgap midpoint", hoverinfo="skip")
    layout = Layout(title=title,
        scene=attr(xaxis=attr(title="k/√n₀", range=[-pi, pi]),
            yaxis=attr(title="x/L", range=[minimum(x), maximum(x)]),
            zaxis=attr(title="ω (rad s⁻¹)", range=[0, 15]),
            camera=attr(eye=attr(x=1.55, y=1.45, z=1.15))),
        width=1100, height=850, margin=attr(l=0, r=0, b=0, t=65),
        legend=attr(x=0.02, y=0.98))
    return PlotlyJS.plot([lower_surface, upper_surface, bandgap_surface], layout)
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
        PlotlyJS.savefig(figure, joinpath(output_directory,
            "frequency_graded_$(name)_3d.html"))
    end
    return output_directory
end