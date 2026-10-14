#=
L1-1 — Bare host: dispersion and Green function.
Pre- and post-processing of the test defined in src/L1_1_bare_host_dispersion.jl
(physics, relation to the paper and observables are documented there).

Study, step 1: one frequency, mesh convergence.
  P1 host (Table 1 of the paper), T_ρ = 0, H = 10 m, ω = 10 rad/s;
  p = 2 and p = 3 with 5, 10, 20, 40, 80 elements per plate wavelength
  (p = 3 at 80 is ≈ 0.8M DOFs and needs more than the 4 GB run limit;
  p = 3 reaches its floor ≈ 1e-6 at 40).

Acceptance (test matrix §0.3, type D):
  - p = 2: observed rate of each error on the two finest meshes ≥ 2(p − 1) − 0.3
    (C/DG plate: energy-norm rate p − 1, smooth functionals 2(p − 1));
  - all orders: error on the finest mesh ≤ 1e-3.

Outputs (relative to the repository root):
  data/sims/L1-1/      one DrWatson record per case (git commit, patch, HydroElasticFEM pin)
  data/exp_pro/L1-1/   errors.csv, checks.csv
  plots/L1-1/          one convergence figure per order (.pdf, .png; dashed:
                       theoretical rate 2(p − 1)) and table (.tex)

Run:
  julia --project -O1 scripts/L1-configA-vertical/L1-1_bare_host_dispersion.jl            # full study
  julia --project -O1 scripts/L1-configA-vertical/L1-1_bare_host_dispersion.jl --smoke    # one coarse case
  julia --project -O1 scripts/L1-configA-vertical/L1-1_bare_host_dispersion.jl --force    # recompute records
=#

using DrWatson
using CairoMakie
using DataFrames
using CSV
using Printf
using TOML
import HydroElasticFEM
using HydroElasticFEM_xxx2026.L1_1

const TEST_ID = "L1-1"
const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# ─── 1. Cases ───────────────────────────────────────────────────────────────

"Cases of the study (`smoke = true`: one coarse case, used by the test suite)."
function cases(; smoke::Bool=false)
    smoke && return [L1_1_Params(elements_per_wavelength=20, order=2, half_window=3.0)]
    return [L1_1_Params(elements_per_wavelength=n, order=p) for p in (2, 3) for n in (5, 10, 20, 40, 80)]
end

# ─── 2. Simulation records (DrWatson) ───────────────────────────────────────

"Record name: the variable parameters, then a hash of the fixed ones."
record_name(p::L1_1_Params) = savename(p; accesses=VARIABLE) * "_fixed=" * fixed_hash(p)

"HydroElasticFEM as pinned in the Manifest (URL, revision, tree hash)."
function hefem_pin()
    manifest = joinpath(REPO, "Manifest.toml")
    isfile(manifest) || return Dict("hefem_source" => "unknown")
    e = only(TOML.parsefile(manifest)["deps"]["HydroElasticFEM"])
    return Dict("hefem_source" => get(e, "repo-url", get(e, "path", "registry")),
                "hefem_rev" => get(e, "repo-rev", "unknown"),
                "hefem_tree" => get(e, "git-tree-sha1", "unknown"),
                "hefem_version" => string(pkgversion(HydroElasticFEM)))
end

"""
Load the record of case `p`, or compute it with `L1_1.run_case` and save it
with all parameters, the git commit and patch of this repository, the
HydroElasticFEM pin, Julia version, wall time and allocations.
"""
function simulate(p::L1_1_Params; root::AbstractString=REPO, force::Bool=false)
    path = joinpath(root, "data", "sims", TEST_ID)
    data, _ = produce_or_load(p, path; filename=record_name, tag=true, gitpath=REPO,
                              storepatch=true, force=force) do p
        stats = @timed run_case(p)
        d = stats.value
        d["params"] = Dict(string(f) => getfield(p, f) for f in fieldnames(L1_1_Params))
        merge!(d, hefem_pin())
        d["julia_version"] = string(VERSION)
        d["wall_time_s"] = stats.time
        d["allocated_bytes"] = stats.bytes
        return d
    end
    return data
end

# ─── 3. Post-processing ─────────────────────────────────────────────────────

"Errors of one case against the references."
function errors(p::L1_1_Params, d)
    γr = d["gamma_r"]
    return (order=p.order, elements_per_wavelength=p.elements_per_wavelength, ω=p.ω,
            ndofs=d["ndofs"], ny=d["ny"], h_over_lambda=d["element_size"] / d["lambda_p"],
            err_k=max(abs(d["k_h_R"] / d["k0"] - 1), abs(d["k_h_L"] / d["k0"] - 1)),
            err_gamma=max(abs(d["gamma_far_R"] - γr), abs(d["gamma_far_L"] - γr)) / γr,
            err_near_re=d["err_near_re"], err_near_im=d["err_near_im"],
            sponge_reflection=d["sponge_reflection"], fit_residual=d["fit_residual"],
            buffer_over_decay=d["buffer"] * d["kappa1"], wall_time_s=d["wall_time_s"])
end

const ERROR_COLUMNS = (:err_k => L"k_h", :err_gamma => L"\gamma_h",
                       :err_near_re => L"\mathrm{Re}\,\eta_h", :err_near_im => L"\mathrm{Im}\,\eta_h")

"Observed order log(eᵢ/eᵢ₊₁)/log(hᵢ/hᵢ₊₁) between the two finest meshes."
last_rate(h, e) = log(e[end - 1] / e[end]) / log(h[end - 1] / h[end])

"Acceptance checks of the study (one row per criterion)."
function checks(df::DataFrame)
    rows = NamedTuple[]
    for g in groupby(sort(df, [:order, :elements_per_wavelength]), :order)
        p = first(g.order)
        for (col, _) in ERROR_COLUMNS
            e = g[end, col]
            push!(rows, (name="$(col), finest mesh, p=$(p)", kind="D", value=e, tolerance=1e-3, passed=e <= 1e-3))
            if p == 2 && nrow(g) >= 2
                rate, theory = last_rate(g.h_over_lambda, g[!, col]), 2(p - 1)
                push!(rows, (name="rate $(col), p=$(p)", kind="D", value=rate, tolerance=theory - 0.3,
                             passed=rate >= theory - 0.3))
            end
        end
    end
    return DataFrame(rows)
end

# ─── 4. Figures and tables ──────────────────────────────────────────────────

"""
One figure per polynomial order: the errors of all quantities against the
mesh size, with the theoretical rate 2(p − 1) (dashed) from the largest
error on the coarsest mesh.
"""
function plot_convergence(df::DataFrame, dir::AbstractString)
    theme = Theme(fontsize=12, palette=(color=[:black, :royalblue, :firebrick, :forestgreen],),
                  Axis=(xgridvisible=false, ygridvisible=false))
    mkpath(dir)
    with_theme(theme) do
        for p in sort(unique(df.order))
            g = sort(filter(r -> r.order == p, df), :elements_per_wavelength)
            n = g.elements_per_wavelength
            fig = Figure(size=(420, 340))
            ax = Axis(fig[1, 1]; xscale=log10, yscale=log10, xlabel="elements per wavelength",
                      ylabel="relative error", title="p = $(p)", xticks=(n, string.(n)))
            for (col, label) in ERROR_COLUMNS
                scatterlines!(ax, n, g[!, col]; label=label)
            end
            rate = 2(p - 1)
            e0 = maximum(g[1, col] for (col, _) in ERROR_COLUMNS)
            lines!(ax, n, e0 .* (n[1] ./ n) .^ rate; color=:gray40, linestyle=:dash,
                   label=L"\mathcal{O}(h^{%$rate})")
            axislegend(ax; position=:lb, framevisible=false)
            save(joinpath(dir, "convergence_single_frequency_p$(p).pdf"), fig)
            save(joinpath(dir, "convergence_single_frequency_p$(p).png"), fig; px_per_unit=300 / 96)
        end
    end
end

function write_table(df::DataFrame, dir::AbstractString)
    mkpath(dir)
    e(x) = @sprintf("%.1e", x)
    open(joinpath(dir, "convergence_single_frequency.tex"), "w") do io
        println(io, "\\begin{table}[htbp]\n\\centering")
        println(io, "\\caption{L1-1: relative errors of the bare plate at \$\\omega = 10\$ rad/s (P1, \$T_\\rho = 0\$, \$H = 10\$ m).}")
        println(io, "\\label{tab:L1-1}\n\\begin{tabular}{rrrcccc}\n\\toprule")
        println(io, "\$p\$ & el./\$\\lambda_p\$ & DOFs & \$k_h\$ & \$\\gamma_h\$ & Re near field & Im near field \\\\\n\\midrule")
        for r in eachrow(sort(df, [:order, :elements_per_wavelength]))
            println(io, join([r.order, r.elements_per_wavelength, r.ndofs, e(r.err_k), e(r.err_gamma),
                              e(r.err_near_re), e(r.err_near_im)], " & "), " \\\\")
        end
        println(io, "\\bottomrule\n\\end{tabular}\n\\end{table}")
    end
end

# ─── 5. Driver ──────────────────────────────────────────────────────────────

"""
    main(; smoke=false, force=false, root=REPO)

Run the study: compute or load every case, write the errors and checks to
`<root>/data/exp_pro/L1-1/`, the figure and table to `<root>/plots/L1-1/`, and
print a summary. Returns `(errors, checks)`.
"""
function main(; smoke::Bool=false, force::Bool=false, root::AbstractString=REPO)
    df = DataFrame([errors(p, simulate(p; root=root, force=force)) for p in cases(; smoke=smoke)])
    ck = checks(df)
    out = joinpath(root, "data", "exp_pro", TEST_ID)
    mkpath(out)
    CSV.write(joinpath(out, "errors.csv"), df)
    CSV.write(joinpath(out, "checks.csv"), ck)
    figdir = joinpath(root, "plots", TEST_ID)
    plot_convergence(df, figdir)
    write_table(df, figdir)

    show(stdout, MIME"text/plain"(), select(df, :order, :elements_per_wavelength, :ndofs, :err_k, :err_gamma,
                                            :err_near_re, :err_near_im); allrows=true)
    println()
    for c in eachrow(ck)
        @printf("%s  %-40s value = %.3g (tolerance %.3g)\n", c.passed ? "PASS" : "FAIL", c.name, c.value, c.tolerance)
    end
    return (errors=df, checks=ck)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(; smoke="--smoke" in ARGS, force="--force" in ARGS)
end
