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

function make_anisotropic_density_panel(title, field; color=:viridis)
    n_wave_numbers = length(field.k_over_sqrt_n0)
    lower_limit = maximum(abs, field.lower_shift)
    upper_limit = maximum(abs, field.upper_shift)
    lower_limits = (-lower_limit, lower_limit)
    upper_limits = (-upper_limit, upper_limit)
    panels = Any[]
    for (branch_name, values, limits) in (("lower shift", field.lower_shift, lower_limits),
        ("upper shift", field.upper_shift, upper_limits))
        for wave_number_index in 1:n_wave_numbers
            heatmap = Plots.heatmap(field.x, field.y,
                values[:, :, wave_number_index];
                xlabel=L"x/L_x", ylabel=L"y/L_y",
                title="$branch_name, k = $(round(field.k_over_sqrt_n0[wave_number_index]; digits=2))",
                color=color, clims=limits, colorbar=true,
                colorbar_title="frequency shift (rad s^-1)", aspect_ratio=1,
                grid=false, fontfamily="Computer Modern", titlefontsize=10,
                guidefontsize=9, tickfontsize=8, colorbar_tickfontsize=8)
            push!(panels, heatmap)
        end
    end
    return Plots.plot(panels...; layout=(2, n_wave_numbers),
        size=(1500, 850), dpi=300, fontfamily="Computer Modern")
end

function make_anisotropic_density_surface(title, table)
    coordinate = table.normalized_positions
    lower_surface = PlotlyJS.surface(x=table.k_over_sqrt_n0, y=coordinate,
        z=table.lower, colorscale=[[0.0, "royalblue"], [1.0, "royalblue"]],
        showscale=false, opacity=0.9, name="lower branch",
        hovertemplate="k/√n₀=%{x:.3f}<br>$(table.axis)/L=%{y:.3f}<br>ω₋=%{z:.3f}<extra></extra>")
    upper_surface = PlotlyJS.surface(x=table.k_over_sqrt_n0, y=coordinate,
        z=table.upper, colorscale=[[0.0, "firebrick"], [1.0, "firebrick"]],
        showscale=false, opacity=0.9, name="upper branch",
        hovertemplate="k/√n₀=%{x:.3f}<br>$(table.axis)/L=%{y:.3f}<br>ω₊=%{z:.3f}<extra></extra>")
    layout = Layout(title=title,
        scene=attr(xaxis=attr(title="k/√n₀", range=[-pi, pi]),
            yaxis=attr(title="$(table.axis)/L", range=[minimum(coordinate), maximum(coordinate)]),
            zaxis=attr(title="ω (rad s⁻¹)", range=[0, 15]),
            camera=attr(eye=attr(x=1.55, y=1.45, z=1.15))),
        width=1100, height=850, margin=attr(l=0, r=0, b=0, t=65),
        legend=attr(x=0.02, y=0.98))
    return PlotlyJS.plot([lower_surface, upper_surface], layout)
end