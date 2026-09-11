function make_anisotropic_density_field_plot(field)
    x_over_Lx = field.x ./ field.Lx
    y_over_Ly = field.y ./ field.Ly
    density_plot = Plots.heatmap(x_over_Lx, y_over_Ly, field.density;
        xlabel=L"x/L_x", ylabel=L"y/L_y", title=L"n(\mathbf{x})/n_0",
        color=:viridis, aspect_ratio=1, colorbar_title=L"n/n_0", grid=false)
    attenuation_plot = Plots.heatmap(x_over_Lx, y_over_Ly, field.attenuation;
        xlabel=L"x/L_x", ylabel=L"y/L_y", title=L"\operatorname{Im} k(\omega_r,\mathbf{x})",
        color=:magma, aspect_ratio=1, colorbar_title=L"\operatorname{Im} k", grid=false)
    return Plots.plot(density_plot, attenuation_plot; layout=(1, 2), size=(1400, 600),
        dpi=300, fontfamily="Computer Modern")
end

function make_anisotropic_optical_depth_plot(aspect_ratios, ratios, profiles)
    ratio_plot = Plots.plot(aspect_ratios, ratios; xlabel=L"L_x/L_y",
        ylabel=L"\mathcal{A}_x/\mathcal{A}_y", label=L"\mathcal{A}_x/\mathcal{A}_y",
        marker=:circle, linewidth=2, color=:black, grid=true)
    Plots.plot!(ratio_plot, aspect_ratios, aspect_ratios; linestyle=:dash,
        color=:gray, label=L"L_x/L_y")
    profile_plot = Plots.plot(profiles.x ./ profiles.Lx, profiles.attenuation_x;
        xlabel=L"x/L_x,\ y/L_y", ylabel=L"\operatorname{Im} k(\omega_r)",
        label=L"x\text{-axis}", linewidth=2, color=:steelblue, grid=true)
    Plots.plot!(profile_plot, profiles.y ./ profiles.Ly, profiles.attenuation_y;
        label=L"y\text{-axis}", linewidth=2, linestyle=:dash, color=:darkorange)
    return Plots.plot(ratio_plot, profile_plot; layout=(1, 2), size=(1400, 600),
        dpi=300, fontfamily="Computer Modern")
end