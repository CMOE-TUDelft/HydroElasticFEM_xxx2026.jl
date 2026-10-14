"""
L1-2 — Flux normalisation (verification test matrix §L1.2).

Cases: P1 host, T_ρ ∈ {0, 9.8}, ω ∈ {2, 5, 10}, p ∈ {2, 3}, mesh-ladder
levels 0–2 (8, 16, 32 elements per λ_p). Observables (see
`src/L1_2_flux_normalisation.jl`): radiative width from the far-field amplitude
and from Im⟨G_b⟩ over the source patch, and the work–flux balance.

Acceptance (type D, §0.3): on the two finest levels the observed rate of each
error is ≥ 2(p − 1) − 0.3, and the far-field γ extrapolated with its
self-convergence rate (three levels, independent of the reference) is within
1e-3 of γ_r = 1/𝒟₀′(k₀). The theoretical rate 2(p − 1) is that of
smooth functionals for the C/DG (C0 interior penalty) plate, whose energy-norm
rate is p − 1.
"""

using .L1_2_FluxNormalisation

const L1_2_ID = "L1-2"

function L1_2_cases(; test_suit::Bool=false)
    test_suit && return [L1_2_params(T_ρ=0.0, ω=10.0, level=0, order=2, half_window=3.0)]
    return vec([L1_2_params(T_ρ=T, ω=ω, level=l, order=p)
                for l in 0:2, p in (2, 3), ω in (2.0, 5.0, 10.0), T in (0.0, 9.8)])
end

function _L1_2_row(c::L1_2_params, d)
    γr = d["gamma_r"]
    return (T_ρ=c.T_ρ, ω=c.ω, order=c.order, level=c.level,
            h_over_lambda=d["element_size"] / d["lambda_p"], gamma_r=γr,
            gamma_far_re=real(d["gamma_far_R"]), gamma_far_im=imag(d["gamma_far_R"]),
            err_far=max(abs(d["gamma_far_R"] - γr), abs(d["gamma_far_L"] - γr)) / γr,
            err_patch=abs(d["gamma_patch"] - γr) / γr,
            err_power=abs(d["P_in"] / d["P_out"] - 1),
            dk_rel=d["k_fit_R"] / d["k0"] - 1,
            sponge_reflection=d["sponge_reflection"], fit_residual=d["fit_residual"])
end

function _L1_2_checks(df::DataFrame)
    checks = AcceptanceCheck[]
    for g in groupby(sort(df, [:T_ρ, :ω, :order, :level]), [:T_ρ, :ω, :order])
        nrow(g) >= 2 || continue
        p = first(g.order)
        tag = @sprintf("T=%g, ω=%g, p=%d", first(g.T_ρ), first(g.ω), p)
        for col in (:err_far, :err_patch, :err_power)
            rate = last(Observables.observed_rates(g.h_over_lambda, g[!, col]))
            add_rate_check!(checks, L1_2_ID, "rate $(col) ($(tag))"; rate=rate, theoretical=2(p - 1))
        end
        nrow(g) >= 3 || continue
        γx = Observables.self_convergence(g.h_over_lambda, complex.(g.gamma_far_re, g.gamma_far_im)).value
        add_check!(checks, L1_2_ID, "extrapolated γ_far ($(tag))";
                   value=real(γx), reference=first(g.gamma_r), tolerance=1e-3, kind="D")
    end
    return checks
end

function _L1_2_plot(df::DataFrame)
    with_theme(VLFS_THEME) do
        fig = Figure(size=(820, 300))
        markers = Dict(2 => :circle, 3 => :rect)
        colors = Dict(2.0 => :black, 5.0 => :royalblue, 10.0 => :firebrick)
        axes = Axis[]
        for (j, T) in enumerate(sort(unique(df.T_ρ)))
            ax = Axis(fig[1, j]; xscale=log10, yscale=log10, xlabel=L"h/\lambda_p",
                      ylabel=j == 1 ? L"|\gamma_h - \gamma_r|/\gamma_r" : "",
                      title=@sprintf("T_ρ = %g m³/s²", T))
            push!(axes, ax)
            for g in groupby(sort(filter(r -> r.T_ρ == T, df), :level), [:ω, :order])
                ω, p = first(g.ω), first(g.order)
                scatterlines!(ax, g.h_over_lambda, g.err_far; color=colors[ω], marker=markers[p],
                              label=@sprintf("ω = %g, p = %d", ω, p))
                lines!(ax, g.h_over_lambda, g.err_patch; color=colors[ω], linestyle=:dot)
            end
        end
        Legend(fig[1, length(axes) + 1], first(axes); framevisible=false, labelsize=9)
        save_figure(fig, L1_2_ID, "gamma_convergence")
    end
end

function _L1_2_table(df::DataFrame)
    finest = combine(groupby(df, [:T_ρ, :ω, :order]), sdf -> sdf[argmax(sdf.level), :])
    tab = select(finest, :T_ρ => "\$T_\\rho\$", :ω => "\$\\omega\$", :order => "\$p\$",
                 :gamma_r => "\$\\gamma_r\$", :gamma_far_re => "\$\\gamma_h\$",
                 :err_far => "far field", :err_patch => "patch", :err_power => "power")
    write_latex_table(L1_2_ID, "gamma_table", tab;
                      caption="L1-2: radiative width from the FE far field, from the source patch and from the work--flux balance (finest mesh).",
                      label="tab:L1-2")
end

function run_L1_2_flux_normalisation(; force=false, make_plots=true, save_csv=true,
                                     verbose=false, vtk_output=false, test_suit=false)
    rows = map(L1_2_cases(; test_suit=test_suit)) do c
        verbose && println("[L1-2] ", savename(c))
        data, _ = produce_or_load_case(run_L1_2, c, L1_2_ID; force=force, verbose=verbose)
        _L1_2_row(c, data)
    end
    df = DataFrame(rows)
    checks = _L1_2_checks(df)
    if save_csv
        mkpath(exp_pro_path(L1_2_ID))
        CSV.write(exp_pro_path(L1_2_ID, "gamma.csv"), df)
        write_checks(L1_2_ID, checks)
    end
    if make_plots
        _L1_2_plot(df)
        _L1_2_table(df)
    end
    return (data=df, checks=checks)
end

register_test(L1_2_ID, run_L1_2_flux_normalisation; cases=L1_2_cases, body=run_L1_2)
