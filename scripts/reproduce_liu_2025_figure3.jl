using HydroElasticFEM_xxx2026
using Plots
using LaTeXStrings

gr()

output_directory = joinpath(@__DIR__, "..", "data", "generated",
    "subsubsec_non_dissipative_homogeneous_isotropic_dispersion")
mkpath(output_directory)

normalized_wave_numbers = collect(range(-pi, pi; length=701))
zero_tension_plate, resonator = liu_2025_parameters()
wave_number_scale = sqrt(resonator.n0)
reference_spacing = inv(wave_number_scale)
pretension = zero_tension_plate.fluid_density * zero_tension_plate.gravity * reference_spacing
pretensioned_plate, _ = liu_2025_parameters(tension=pretension)

cases = [
    ("bare_T0", L"\mathrm{Bare\ plate},\ T=0", zero_tension_plate, false),
    ("LRH_T0", L"\mathrm{LRH\ plate},\ T=0", zero_tension_plate, true),
    ("bare_Tpositive", L"\mathrm{Bare\ plate},\ T>0", pretensioned_plate, false),
    ("LRH_Tpositive", L"\mathrm{LRH\ plate},\ T>0", pretensioned_plate, true),
]

function physical_wavenumber(normalized_k)
    return max(abs(normalized_k) * wave_number_scale, eps(Float64))
end

function case_curves(plate, include_resonator)
    physical_wave_numbers = physical_wavenumber.(normalized_wave_numbers)
    bare = [bare_dispersion(plate, k) for k in physical_wave_numbers]
    if include_resonator
        lower = [metaplate_dispersion(plate, resonator, k).lower for k in physical_wave_numbers]
        upper = [metaplate_dispersion(plate, resonator, k).upper for k in physical_wave_numbers]
        return bare, lower, upper
    end
    return bare, nothing, nothing
end

function add_bandgap_shading!(panel, lower, upper)
    gap_low = maximum(lower)
    gap_high = minimum(upper)
    if gap_low < gap_high
        gap_shape = Shape([-pi, pi, pi, -pi],
            [gap_low, gap_low, gap_high, gap_high])
        plot!(panel, gap_shape; label=L"\mathrm{band\ gap}", color=:gray, fillalpha=0.35,
            linecolor=:transparent)
    end
end

function make_panel(title, plate, include_resonator; y_limits=nothing)
    bare, lower, upper = case_curves(plate, include_resonator)
    resolved_y_limits = isnothing(y_limits) ? (0.0, 15.0) : y_limits
    panel = plot(normalized_wave_numbers, bare;
        label=include_resonator ? L"\mathrm{Bare\ reference}" : L"\mathrm{Bare\ plate}",
        color=:black, linewidth=2, xlabel=L"k/\sqrt{n_0}", ylabel=L"\omega\ (\mathrm{rad\ s^{-1}})",
        title=title, legend=:top, grid=true, xlims=(-pi, pi), ylims=resolved_y_limits)
    if include_resonator
        add_bandgap_shading!(panel, lower, upper)
        plot!(panel, normalized_wave_numbers, lower;
            label=L"\mathrm{LRH\ plate\ lower}", color=:royalblue, linewidth=2)
        plot!(panel, normalized_wave_numbers, upper;
            label=L"\mathrm{LRH\ plate\ upper}", color=:firebrick, linewidth=2)
        plot!(panel, [-pi, pi],
            [resonator.natural_frequency, resonator.natural_frequency];
            label=L"\omega_r", color=:gray, linestyle=:dash)
    end
    return panel
end

function bandgap_limits(plate)
    _, lower, upper = case_curves(plate, true)
    gap_low = maximum(lower)
    gap_high = minimum(upper)
    margin = 0.10 * (gap_high - gap_low)
    return (max(0.0, gap_low - margin), min(15.0, gap_high + margin))
end

for (case_name, title, plate, include_resonator) in cases
    figure = make_panel(title, plate, include_resonator)
    savefig(figure, joinpath(output_directory, "$(case_name).png"))
end

zoomed_plots = [
    make_panel(L"\mathrm{LRH\ plate},\ T=0", zero_tension_plate, true;
        y_limits=bandgap_limits(zero_tension_plate)),
    make_panel(L"\mathrm{LRH\ plate},\ T>0", pretensioned_plate, true;
        y_limits=bandgap_limits(pretensioned_plate)),
]
for (case_name, zoomed_figure) in zip(("LRH_T0", "LRH_Tpositive"), zoomed_plots)
    savefig(zoomed_figure,
    joinpath(output_directory, "$(case_name)_bandgap_zoom.png"))
end

open(joinpath(output_directory, "dispersion.csv"), "w") do io
    println(io, "k_over_sqrt_n0,bare_T0,lrh_lower_T0,lrh_upper_T0,bare_Tpositive,lrh_lower_Tpositive,lrh_upper_Tpositive")
    for (normalized_k, physical_k) in zip(normalized_wave_numbers,
        physical_wavenumber.(normalized_wave_numbers))
        no_tension = metaplate_dispersion(zero_tension_plate, resonator, physical_k)
        with_tension = metaplate_dispersion(pretensioned_plate, resonator, physical_k)
        println(io, join((normalized_k, bare_dispersion(zero_tension_plate, physical_k), no_tension.lower,
            no_tension.upper, bare_dispersion(pretensioned_plate, physical_k),
            with_tension.lower, with_tension.upper), ","))
    end
end

println("Wrote four standalone case figures in ", output_directory)
println("Wrote two LRH band-gap close-ups in ", output_directory)
println("Wrote ", joinpath(output_directory, "dispersion.csv"))