function make_frequency_graded_panel(
    title,
    table;
    colors=palette(:viridis, length(table.x_over_L)),
    y_limits=(0.0, 15.0),
)
    panel = plot(; xlabel=L"k/\sqrt{n_0}", ylabel=L"\omega\ (\mathrm{rad\ s^{-1}})",
        title=title, legend=:top, grid=true, xlims=(-pi, pi), ylims=y_limits,
        size=(900, 650), dpi=300, fontfamily="Computer Modern")
    for (position_index, position) in enumerate(table.x_over_L)
        label = L"x/L=%$(round(position; digits=2))"
        plot!(panel, table.k_over_sqrt_n0, table.lower[position_index, :];
            color=colors[position_index], linewidth=1.5, label=label)
        plot!(panel, table.k_over_sqrt_n0, table.upper[position_index, :];
            color=colors[position_index], linewidth=1.5, linestyle=:dash, label=false)
    end
    return panel
end
