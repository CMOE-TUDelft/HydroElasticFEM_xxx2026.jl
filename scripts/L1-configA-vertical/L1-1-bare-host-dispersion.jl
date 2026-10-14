"""
L1-1 — Bare-host dispersion and Green function (verification test matrix §L1.1).

Step 1 (this file): convergence at one frequency. P1 host, T_ρ = 0,
ω = 10 rad/s, H = 10 m; 20, 40, 80 elements per plate wavelength for p = 2
and 20, 40 for p = 3 (finer meshes exceed ≈ 0.8M DOFs, beyond the 4 GB run
limit with the κ₁-based buffer; p = 3 reaches its floor by 40). The sponges
are ≥ 5/κ₁ from the source (builder default).
Observables (`src/L1_1_bare_host_dispersion.jl`): relative errors of the
discrete wavenumber k_h, of the far-field radiative width γ, and of the near
field (real and imaginary parts) against the patch Green function.

Acceptance (type D, §0.3): for p = 2, the observed rate of each error on the
two finest meshes is ≥ 2(p − 1) − 0.3 (C/DG plate, smooth functionals); for
both orders, the error on the finest mesh is ≤ 1e-3.
"""

using .L1_1_BareHostDispersion

const L1_1_ID = "L1-1"

function L1_1_cases(; test_suit::Bool=false)
    test_suit && return [L1_1_params(elements_per_wavelength=20, half_window=3.0)]
    return [[L1_1_params(elements_per_wavelength=n, order=2) for n in (20, 40, 80)];
            [L1_1_params(elements_per_wavelength=n, order=3) for n in (20, 40)]]
end

function _L1_1_row(c::L1_1_params, d)
    γr = d["gamma_r"]
    return (T_ρ=c.T_ρ, H=c.H, ω=c.ω, order=c.order, elements_per_wavelength=c.elements_per_wavelength,
            h_over_lambda=d["element_size"] / d["lambda_p"], ny=d["ny"], ndofs=d["ndofs"],
            err_k=max(abs(d["k_h_R"] / d["k0"] - 1), abs(d["k_h_L"] / d["k0"] - 1)),
            err_gamma=max(abs(d["gamma_far_R"] - γr), abs(d["gamma_far_L"] - γr)) / γr,
            err_near_re=d["err_near_re"], err_near_im=d["err_near_im"],
            err_near_re_evan=d["err_near_re_evan"], err_near_im_evan=d["err_near_im_evan"],
            sponge_distance=d["sponge_distance"], kappa1=d["kappa1"],
            sponge_reflection=d["sponge_reflection"], fit_residual=d["fit_residual"],
            wall_time_s=d["wall_time_s"])
end

const _L1_1_ERRORS = (:err_k => L"k_h", :err_gamma => L"\gamma_h",
                      :err_near_re => L"\mathrm{Re}\,\eta_h", :err_near_im => L"\mathrm{Im}\,\eta_h")

function _L1_1_checks(df::DataFrame)
    checks = AcceptanceCheck[]
    for g in groupby(sort(df, [:T_ρ, :H, :ω, :order, :elements_per_wavelength]), [:T_ρ, :H, :ω, :order])
        p = first(g.order)
        tag = @sprintf("T=%g, H=%g, ω=%g, p=%d", first(g.T_ρ), first(g.H), first(g.ω), p)
        for (col, _) in _L1_1_ERRORS
            add_check!(checks, L1_1_ID, "$(col) finest ($(tag))";
                       value=g[end, col], reference=0.0, tolerance=1e-3, kind="D")
            if p == 2 && nrow(g) >= 2
                rate = last(Observables.observed_rates(g.h_over_lambda, g[!, col]))
                add_rate_check!(checks, L1_1_ID, "rate $(col) ($(tag))"; rate=rate, theoretical=2(p - 1))
            end
        end
    end
    return checks
end

function _L1_1_plot(df::DataFrame)
    with_theme(VLFS_THEME) do
        fig = Figure(size=(820, 320))
        orders = sort(unique(df.order))
        for (j, p) in enumerate(orders)
            g = sort(filter(r -> r.order == p, df), :elements_per_wavelength)
            ax = Axis(fig[1, j]; xscale=log10, yscale=log10, xlabel="elements per wavelength",
                      ylabel=j == 1 ? "relative error" : "", title="p = $(p)",
                      xticks=(g.elements_per_wavelength, string.(g.elements_per_wavelength)))
            for (col, label) in _L1_1_ERRORS
                scatterlines!(ax, g.elements_per_wavelength, g[!, col]; label=label)
            end
            j == length(orders) && axislegend(ax; position=:lb, framevisible=false, labelsize=10)
        end
        save_figure(fig, L1_1_ID, "convergence_single_frequency")
    end
end

function _L1_1_table(df::DataFrame)
    tab = select(sort(df, [:order, :elements_per_wavelength]), :order => "\$p\$",
                 :elements_per_wavelength => "el./\$\\lambda_p\$", :ndofs => "DOFs",
                 :err_k => "\$k_h\$", :err_gamma => "\$\\gamma_h\$",
                 :err_near_re => "Re near field", :err_near_im => "Im near field")
    write_latex_table(L1_1_ID, "convergence_single_frequency", tab;
                      caption="L1-1: relative errors of the bare plate at \$\\omega = 10\$ rad/s (P1, \$T_\\rho = 0\$, \$H = 10\$ m).",
                      label="tab:L1-1")
end

function run_L1_1_bare_host_dispersion(; force=false, make_plots=true, save_csv=true,
                                       verbose=false, vtk_output=false, test_suit=false)
    rows = map(L1_1_cases(; test_suit=test_suit)) do c
        verbose && println("[L1-1] ", savename(c))
        data, _ = produce_or_load_case(run_L1_1, c, L1_1_ID; force=force, verbose=verbose)
        _L1_1_row(c, data)
    end
    df = DataFrame(rows)
    checks = _L1_1_checks(df)
    if save_csv
        mkpath(exp_pro_path(L1_1_ID))
        CSV.write(exp_pro_path(L1_1_ID, "convergence_single_frequency.csv"), df)
        write_checks(L1_1_ID, checks)
    end
    if make_plots
        _L1_1_plot(df)
        _L1_1_table(df)
    end
    return (data=df, checks=checks)
end

register_test(L1_1_ID, run_L1_1_bare_host_dispersion; cases=L1_1_cases, body=run_L1_1)
