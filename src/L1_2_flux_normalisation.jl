"""
L1-2 — Flux normalisation (`eq:lat_flux`, `eq:lat_gamma_r`).

A density-normalised force F on the bare P1 plate (configuration A) radiates
waves of amplitude A = iF sinc(k₀b/2)/𝒟₀′(k₀) to both sides (patch of width b,
`References.green_1d_patch`). Three FE measures of the radiative width
γ_r = Im G(0) = 1/𝒟₀′(k₀) are returned:

- far field: the source amplitude from the waves fitted on probe windows on
  either side, γ = A/(iF sinc(k₀b/2)), with the waves entering from the
  sponges removed (the source is transparent) and the discrete wavenumber
  fitted, so that the dispersion error (L1-1) does not enter the phase of A;
- near field: Im⟨G_b⟩ over the patch, since Im G(x) = cos(k₀x)/𝒟₀′ for a lossless
  host gives Im⟨G_b⟩ = sinc²(k₀b/2)/𝒟₀′;
- power: work done by the force, ½ρω Im(F̄⟨η⟩), against the net flux
  ρω/4 𝒟₀′ Σ(|A_out|² − |A_in|²) through both probe windows.
"""
module L1_2_FluxNormalisation

using Parameters
using ..ParameterSets: P1
using ..References: k0, dD0, gamma_r, host_roots
using ..FEBuilders: solve_configA, ladder_mesh
using ..Observables: fit_waves, patch_average

export L1_2_params, run_L1_2

"""
Case of L1-2. `level` is the mesh-ladder level (`FEBuilders.ladder_mesh`);
`half_window` and `probe_offset` are in plate wavelengths: the probes span
[x_s + probe_offset λ_p, x_s + half_window λ_p] on either side of the source.
"""
@with_kw struct L1_2_params
    T_ρ::Float64 = 0.0
    ω::Float64 = 10.0
    level::Int = 1
    order::Int = 2
    half_window::Float64 = 10.0
    probe_offset::Float64 = 2.0
    n_modes::Int = 6
end

"Run one L1-2 case; returns the raw FE measures and their references."
function run_L1_2(p::L1_2_params)
    @unpack T_ρ, ω, level, order, half_window, probe_offset, n_modes = p
    h = P1(T_ρ=T_ρ).host
    F = 1.0
    kp = k0(h, ω)
    λ_p = 2π / kp
    mesh = ladder_mesh(h, ω, level; order=order)
    sol = solve_configA(h, ω; window_length=2half_window * λ_p, source_x=half_window * λ_p,
                        F=F, mesh=mesh)
    lay = sol.layout
    xs, b = lay.source_x, lay.source_width
    roots = host_roots(h, ω)

    right = fit_waves(sol.ηh, roots, xs + probe_offset * λ_p, lay.window[2]; x_ref=xs, n_modes=n_modes)
    left = fit_waves(sol.ηh, roots, lay.window[1], xs - probe_offset * λ_p; x_ref=xs, n_modes=n_modes)
    out_R, in_R = right.Aplus, right.Aminus
    out_L, in_L = left.Aminus, left.Aplus

    s = sin(kp * b / 2) / (kp * b / 2)
    γ_R = (out_R - in_L) / (im * s * F)
    γ_L = (out_L - in_R) / (im * s * F)

    ηavg = patch_average(sol.ηh, xs - b / 2, xs + b / 2, round(Int, b / lay.element_size); npts=order + 2)
    γ_patch = imag((ηavg - s * (in_L + in_R)) / F) / s^2

    d = dD0(h, kp, ω)
    P_in = h.ρ * ω / 2 * imag(conj(F) * ηavg)
    P_out = h.ρ * ω / 4 * d * (abs2(out_R) - abs2(in_R) + abs2(out_L) - abs2(in_L))

    return Dict{String,Any}(
        "gamma_r" => gamma_r(h, ω; dim=2),
        "gamma_far_R" => γ_R, "gamma_far_L" => γ_L, "gamma_patch" => γ_patch,
        "P_in" => P_in, "P_out" => P_out,
        "sponge_reflection" => max(abs(in_R / out_R), abs(in_L / out_L)),
        "fit_residual" => max(right.residual, left.residual),
        "k_fit_R" => right.k, "k_fit_L" => left.k,
        "k0" => kp, "lambda_p" => λ_p, "b" => b, "element_size" => lay.element_size,
        "ny" => mesh.ny, "vertical_ratio" => mesh.vertical_ratio,
    )
end

end # module L1_2_FluxNormalisation
