function run_non_dissipative_homogeneous_isotropic_lrh_plate(; output_directory=joinpath(
        @__DIR__, "..", "..", "..", "data", "generated",
        "dispersion_theory", "homogeneous_isotropic"), n_points=701)
    gr()
    mkpath(output_directory)
    normalized_wave_numbers = collect(range(-pi, pi; length=n_points))
    zero_tension_plate, resonator = liu_2025_parameters()
    reference_spacing = inv(sqrt(resonator.n0))
    pretension = zero_tension_plate.fluid_density * zero_tension_plate.gravity * reference_spacing
    pretensioned_plate, _ = liu_2025_parameters(tension=pretension)
    cases = (("bare_T0", L"\mathrm{Bare\ plate},\ T=0", zero_tension_plate, false),
        ("LRH_T0", L"\mathrm{LRH\ plate},\ T=0", zero_tension_plate, true),
        ("bare_Tpositive", L"\mathrm{Bare\ plate},\ T>0", pretensioned_plate, false),
        ("LRH_Tpositive", L"\mathrm{LRH\ plate},\ T>0", pretensioned_plate, true))
    tables = (zero_tension=dispersion_table(zero_tension_plate, resonator,
        normalized_wave_numbers), pretensioned=dispersion_table(pretensioned_plate,
        resonator, normalized_wave_numbers))
    for (case_name, title, plate, include_resonator) in cases
        table = dispersion_table(plate, resonator, normalized_wave_numbers;
            include_resonator=include_resonator)
        savefig(make_non_dissipative_homogeneous_isotropic_lrh_plate_plot(title, table,
            resonator), joinpath(output_directory, "$(case_name).png"))
    end
    savefig(make_non_dissipative_homogeneous_isotropic_lrh_plate_plot(
        L"\mathrm{LRH\ plate},\ T=0", tables.zero_tension, resonator;
        y_limits=(9.5,10.5)), joinpath(output_directory, "LRH_T0_bandgap_zoom.png"))
    savefig(make_non_dissipative_homogeneous_isotropic_lrh_plate_plot(
        L"\mathrm{LRH\ plate},\ T>0", tables.pretensioned, resonator;
        y_limits=(9.5,10.5)),
        joinpath(output_directory, "LRH_Tpositive_bandgap_zoom.png"))
    write_dispersion_csv(joinpath(output_directory, "dispersion.csv"), tables)
    return output_directory
end