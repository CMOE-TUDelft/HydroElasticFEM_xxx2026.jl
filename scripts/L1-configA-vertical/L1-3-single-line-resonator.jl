"""
L1-3 — Single line resonator (verification test matrix §L1.3).

Cases (P1 host, m_ρ,r = 0.01, ω_r = 10):
- convergence ladder: T_ρ ∈ {0, 9.8}, c ∈ {0, c_crit}, resonator on a node;
  c = 0 also at an element midpoint. Levels 0–2 for p = 3 and 1–3 for p = 2:
  P2 at level 0 (8 elements per λ_p) does not resolve the resonance (the
  renormalised detuning is not even monotone in ω), so the P2 ladder starts
  at 16 elements per λ_p;
- response curves: T_ρ ∈ {0, 9.8}, c/c_crit ∈ {0, ½, 1, 2}, level 1, p = 3,
  101 frequencies in ω_r(1 ± 5/Q_rad).

Observables and references (`src/L1_3_single_line_resonator.jl`): renormalised
frequency ω̃ (`eq:lat_renorm_single`), R = −1 and T = 0 at ω̃ for c = 0,
𝒜 = ½ at critical coupling (`eq:lat_absorbed_1d`), half-power width, and the
resonator response q/η.

Acceptance (type D, §0.3): ω̃ and half-power width extrapolated with their
self-convergence rate (three levels, independent of the reference) within
1e-3 of the reference, |R + 1| and |T| at ω̃ (c = 0) and |𝒜 − ½|
(c = c_crit) at most 1e-3 on the finest p = 3 mesh, and rate ≥ 2(p − 1) − 0.3
for the ω̃ error. The width is also compared with 2/Q_rad (type A, reported):
Δω/ω̃ ≈ (1 + c/c_crit)/Q_rad, i.e. 2/Q_rad at critical coupling.
"""

using .L1_3_SingleLineResonator

const L1_3_ID = "L1-3"

function L1_3_cases(; test_suit::Bool=false)
    test_suit && return [L1_3_params(level=0, order=3, n_omega=3)]
    levels(p) = p == 2 ? (1:3) : (0:2)
    ladder = [L1_3_params(T_ρ=T, c_factor=c, level=l, order=p)
              for T in (0.0, 9.8) for c in (0.0, 1.0) for p in (2, 3) for l in levels(p)]
    midelement = [L1_3_params(T_ρ=T, snap=:midelement, level=l, order=p)
                  for T in (0.0, 9.8) for p in (2, 3) for l in levels(p)]
    curves = [L1_3_params(T_ρ=T, c_factor=c, level=1, order=3, n_omega=101)
              for c in (0.0, 0.5, 1.0, 2.0), T in (0.0, 9.8)]
    return [ladder; midelement; vec(curves)]
end

"""
Reference resonance of the case: root of Im(1/R) = −x/γ_r of the exact
single-resonator solution with the case's dashpot (the record's
`omega_peak_ref` is the undamped renormalised frequency, which only defines
c_crit).
"""
function _L1_3_reference_peak(c::L1_3_params, ω̃0::Real, Q::Real)
    h, r, _, _ = L1_3_SingleLineResonator.L1_3_resonator(c)
    f(ω) = imag(1 / References.single_resonator_1d(h, r, ω).R)
    Δ = 4c.ω_r / Q
    return References.find_zero(f, (ω̃0 - Δ, ω̃0 + Δ))
end

function _L1_3_row(c::L1_3_params, d)
    width_h = -(d["omega_half_h"]...)
    width_ref = -(d["omega_half_ref"]...)
    ω̃_ref = _L1_3_reference_peak(c, d["omega_peak_ref"], d["Q_rad"])
    return (T_ρ=c.T_ρ, c_factor=c.c_factor, snap=String(c.snap), order=c.order, level=c.level,
            n_omega=c.n_omega, h_over_lambda=d["element_size"] / d["lambda_p"],
            omega_peak_h=d["omega_peak_h"], omega_peak_ref=ω̃_ref,
            err_omega_peak=abs(d["omega_peak_h"] / ω̃_ref - 1),
            width_h=abs(width_h), width_ref=abs(width_ref),
            err_width=abs(width_h / width_ref - 1),
            # Γ/γ_r ≈ 1 + c/c_crit, so Δω/ω̃ ≈ (1 + c/c_crit)/Q_rad (2/Q_rad at critical coupling)
            width_Qrad=(1 + c.c_factor) / d["Q_rad"] * ω̃_ref,
            R_peak_re=real(d["R_peak_h"]), R_peak_im=imag(d["R_peak_h"]),
            abs_T_peak=abs(d["T_peak_h"]), A_peak_h=d["A_peak_h"], A_peak_ref=d["A_peak_ref"],
            err_q=abs(d["q_over_eta_peak_h"] / d["q_over_eta_peak_ref"] - 1),
            Q_rad=d["Q_rad"], offset_in_element=d["offset_in_element"],
            fit_residual=d["fit_residual_peak"], dk_rel=d["dispersion_error"])
end

function _L1_3_checks(df::DataFrame)
    checks = AcceptanceCheck[]
    ladder = filter(r -> r.n_omega == 0, df)
    for g in groupby(sort(ladder, [:T_ρ, :c_factor, :snap, :order, :level]), [:T_ρ, :c_factor, :snap, :order])
        nrow(g) >= 2 || continue
        p = first(g.order)
        tag = @sprintf("T=%g, c/c_crit=%g, %s, p=%d", first(g.T_ρ), first(g.c_factor), first(g.snap), p)
        hs = g.h_over_lambda
        rate = last(Observables.observed_rates(hs, g.err_omega_peak))
        add_rate_check!(checks, L1_3_ID, "rate ω̃ ($(tag))"; rate=rate, theoretical=2(p - 1))
        for (col, ref) in ((:omega_peak_h, :omega_peak_ref), (:width_h, :width_ref))
            nrow(g) >= 3 || continue
            x = Observables.self_convergence(hs, g[!, col]).value
            add_check!(checks, L1_3_ID, "extrapolated $(col) ($(tag))";
                       value=x, reference=g[end, ref], tolerance=1e-3, kind="D")
        end
        p == 3 || continue
        fin = g[end, :]
        if fin.c_factor == 0
            add_check!(checks, L1_3_ID, "R(ω̃) = −1 ($(tag))";
                       value=hypot(fin.R_peak_re + 1, fin.R_peak_im), reference=0.0, tolerance=1e-3, kind="D")
            add_check!(checks, L1_3_ID, "T(ω̃) = 0 ($(tag))";
                       value=fin.abs_T_peak, reference=0.0, tolerance=1e-3, kind="D")
        elseif fin.c_factor == 1
            add_check!(checks, L1_3_ID, "𝒜(ω̃) = ½ at critical coupling ($(tag))";
                       value=fin.A_peak_h, reference=0.5, tolerance=1e-3, kind="D")
        end
        add_check!(checks, L1_3_ID, "q/η at ω̃ ($(tag))";
                   value=fin.err_q, reference=0.0, tolerance=1e-3, kind="D")
    end
    return checks
end

function _L1_3_curves(records)
    rows = NamedTuple[]
    for (c, d) in records
        c.n_omega > 0 || continue
        for i in eachindex(d["omega"])
            push!(rows, (T_ρ=c.T_ρ, c_factor=c.c_factor, omega=d["omega"][i],
                         abs_R_h=abs(d["R_h"][i]), abs_T_h=abs(d["T_h"][i]), A_h=d["A_h"][i],
                         abs_R_ref=abs(d["R_ref"][i]), abs_T_ref=abs(d["T_ref"][i]), A_ref=d["A_ref"][i]))
        end
    end
    return DataFrame(rows)
end

function _L1_3_plot(curves::DataFrame)
    nrow(curves) == 0 && return
    with_theme(VLFS_THEME) do
        Ts = sort(unique(curves.T_ρ))
        fig = Figure(size=(820, 520))
        colors = Dict(0.0 => :black, 0.5 => :royalblue, 1.0 => :firebrick, 2.0 => :forestgreen)
        for (j, T) in enumerate(Ts)
            axR = Axis(fig[1, j]; ylabel=j == 1 ? L"|R|" : "", title=@sprintf("T_ρ = %g m³/s²", T))
            axA = Axis(fig[2, j]; xlabel=L"\omega\ \mathrm{(rad/s)}", ylabel=j == 1 ? L"\mathcal{A}" : "")
            hlines!(axA, [0.5]; color=:gray, linestyle=:dash)
            for g in groupby(sort(filter(r -> r.T_ρ == T, curves), :omega), :c_factor)
                cf = first(g.c_factor)
                col = colors[cf]
                lines!(axR, g.omega, g.abs_R_ref; color=col, label=@sprintf("c/c_crit = %g", cf))
                scatter!(axR, g.omega[1:4:end], g.abs_R_h[1:4:end]; color=col, markersize=5)
                lines!(axA, g.omega, g.A_ref; color=col)
                scatter!(axA, g.omega[1:4:end], g.A_h[1:4:end]; color=col, markersize=5)
            end
            j == 1 && axislegend(axR; position=:rt, framevisible=false, labelsize=9)
        end
        save_figure(fig, L1_3_ID, "single_resonator_response")
    end
end

function _L1_3_table(df::DataFrame)
    fin = filter(r -> r.order == maximum(df.order) && r.snap == "node", df)
    fin = combine(groupby(fin, [:T_ρ, :c_factor]), sdf -> sdf[argmax(sdf.level), :])
    tab = select(fin, :T_ρ => "\$T_\\rho\$", :c_factor => "\$c/c_{\\rm crit}\$",
                 :omega_peak_h => "\$\\tilde\\omega_{r,h}\$", :omega_peak_ref => "\$\\tilde\\omega_r\$",
                 :width_h => "\$\\Delta\\omega_h\$", :width_ref => "\$\\Delta\\omega\$",
                 :width_Qrad => "\$\\Delta\\omega(Q_{\\rm rad})\$", :A_peak_h => "\$\\mathcal{A}_h(\\tilde\\omega)\$",
                 :Q_rad => "\$Q_{\\rm rad}\$")
    write_latex_table(L1_3_ID, "single_resonator_table", tab;
                      caption="L1-3: renormalised resonance, half-power width and absorbed fraction of one line resonator (FE, \$p = 3\$, finest mesh) against the exact single-resonator solution.",
                      label="tab:L1-3")
end

function run_L1_3_single_line_resonator(; force=false, make_plots=true, save_csv=true,
                                        verbose=false, vtk_output=false, test_suit=false)
    records = map(L1_3_cases(; test_suit=test_suit)) do c
        verbose && println("[L1-3] ", savename(c))
        data, _ = produce_or_load_case(run_L1_3, c, L1_3_ID; force=force, verbose=verbose)
        (c, data)
    end
    df = DataFrame([_L1_3_row(c, d) for (c, d) in records])
    curves = _L1_3_curves(records)
    checks = _L1_3_checks(df)
    if save_csv
        mkpath(exp_pro_path(L1_3_ID))
        CSV.write(exp_pro_path(L1_3_ID, "peak_and_width.csv"), df)
        nrow(curves) > 0 && CSV.write(exp_pro_path(L1_3_ID, "response_curves.csv"), curves)
        write_checks(L1_3_ID, checks)
    end
    if make_plots
        _L1_3_plot(curves)
        _L1_3_table(df)
    end
    return (data=df, curves=curves, checks=checks)
end

register_test(L1_3_ID, run_L1_3_single_line_resonator; cases=L1_3_cases, body=run_L1_3)
