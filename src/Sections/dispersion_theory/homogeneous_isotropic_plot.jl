function add_bandgap_shading!(panel, table)
    gap_low = maximum(table.lower)
    gap_high = minimum(table.upper)
    if gap_low < gap_high
        k_min, k_max = extrema(table.k_over_sqrt_n0)
        gap_shape = Shape([k_min, k_max, k_max, k_min],
            [gap_low, gap_low, gap_high, gap_high])
        Plots.plot!(panel, gap_shape; label=L"\mathrm{band\ gap}", color=:gray,
            fillalpha=0.35, linecolor=:transparent)
    end
end

function make_non_dissipative_homogeneous_isotropic_lrh_plate_plot(
    title, table, resonator; y_limits=(0.0, 15.0))
    panel = Plots.plot(table.k_over_sqrt_n0, table.bare;
        label=table.lower === nothing ? L"\mathrm{Bare\ plate}" : L"\mathrm{Bare\ reference}",
        color=:black, linewidth=2, xlabel=L"k/\sqrt{n_0}",
        ylabel=L"\omega\ (\mathrm{rad\ s^{-1}})", title=title,
        legend=:top, grid=true, xlims=extrema(table.k_over_sqrt_n0), ylims=y_limits)
    if table.lower !== nothing
        add_bandgap_shading!(panel, table)
        Plots.plot!(panel, table.k_over_sqrt_n0, table.lower;
            label=L"\mathrm{LRH\ plate\ lower}", color=:royalblue, linewidth=2)
        Plots.plot!(panel, table.k_over_sqrt_n0, table.upper;
            label=L"\mathrm{LRH\ plate\ upper}", color=:firebrick, linewidth=2)
        k_min, k_max = extrema(table.k_over_sqrt_n0)
        Plots.plot!(panel, [k_min, k_max], [resonator.natural_frequency,
            resonator.natural_frequency]; label=L"\omega_r", color=:gray, linestyle=:dash)
    end
    return panel
end