function run_anisotropic_density_dispersion(; output_directory=joinpath(
        @__DIR__, "..", "..", "..", "data", "generated", "dispersion_theory",
    "anisotropic_density"), aspect_ratios=[1.0, 2.0, 4.0, 6.0],
    damping_ratio=0.05, n_points=701, n_positions=31)
    Plots.gr()
    mkpath(output_directory)
    reference_plate, resonator = liu_2025_parameters()
    pretension = reference_plate.fluid_density * reference_plate.gravity /
        sqrt(resonator.n0)
    plate = PlateParameters(reference_plate.fluid_density, reference_plate.gravity,
        reference_plate.depth, reference_plate.mass_density, pretension, 0.0)
    field = anisotropic_density_field(plate, resonator; Lx=4.0, Ly=1.0,
        damping_ratio=damping_ratio)
    field_figure = make_anisotropic_density_field_plot(field)
    Plots.savefig(field_figure, joinpath(output_directory, "anisotropic_density_field.png"))

    optical_depths = [anisotropic_optical_depths(plate, resonator;
        Lx=aspect_ratio, Ly=1.0, damping_ratio=damping_ratio) for aspect_ratio in aspect_ratios]
    ratios = [depth.A_x / depth.A_y for depth in optical_depths]
    profiles = optical_depths[findfirst(==(4.0), aspect_ratios)]
    optical_figure = make_anisotropic_optical_depth_plot(aspect_ratios, ratios, profiles)
    Plots.savefig(optical_figure,
        joinpath(output_directory, "anisotropic_density_optical_depth.png"))

    normalized_wave_numbers = collect(range(-pi, pi; length=n_points))
    normalized_positions = collect(range(-4.0, 4.0; length=n_positions))
    tables = (x=anisotropic_density_table(plate, resonator, normalized_wave_numbers,
            normalized_positions; Lx=4.0, Ly=1.0, axis=:x),
        y=anisotropic_density_table(plate, resonator, normalized_wave_numbers,
            normalized_positions; Lx=4.0, Ly=1.0, axis=:y))
    field_wave_numbers = [0.25, 0.75, 1.5, 2.5]
    dispersion_field = anisotropic_dispersion_field(plate, resonator,
        field_wave_numbers; Lx=4.0, Ly=1.0)
    dispersion_figure = make_anisotropic_density_panel(
        L"Anisotropic density: local frequency shift, T>0", dispersion_field;
        color=:RdBu)
    Plots.savefig(dispersion_figure,
        joinpath(output_directory, "anisotropic_density_dispersion.png"))
    for axis in (:x, :y)
        figure = make_anisotropic_density_surface(
            "Anisotropic density, $(axis)-slice, T > 0", tables[axis])
        PlotlyJS.savefig(figure, joinpath(output_directory,
            "anisotropic_density_$(axis)_3d.html"))
    end
    return output_directory
end