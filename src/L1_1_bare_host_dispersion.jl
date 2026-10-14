"""
L1-1 — Bare-host dispersion and Green function (`eq:lat_symbol`,
`eq:lat_green_1d`, `eq:lat_flux`).

The bare P1 plate (configuration A) is driven by a density-normalised force F
spread over a patch of width b at the centre of the window. At one frequency
and one mesh the run returns:

- the discrete wavenumber k_h, fitted on probe windows on either side of the
  source (`Observables.fit_waves`), against k₀;
- the radiative width from the far-field amplitude, γ = A/(iF sinc(k₀b/2)),
  with the waves entering from the sponges removed, against γ_r = 1/𝒟₀′(k₀);
- the near field between the patch edge and the probe windows, after removing
  the propagating waves entering from the sponges, against F·G_b(x − x_s)
  (`References.green_1d_patch`), with the real and imaginary parts reported
  separately (Re G is the part of the Green function that sets the
  renormalised resonance of L1-3). The error after also removing the fitted
  incoming evanescent modes (`Observables.incoming_field`) is returned as a
  diagnostic (`*_evan`): on 3 λ_p probe windows the slowest evanescent mode
  entering from the end and the one leaving the source are nearly collinear,
  so that correction is ill-conditioned (≈ 1e-2 errors at ω = 10); the
  κ₁-based buffer removes the evanescent reflection instead.
"""
module L1_1_BareHostDispersion

using Parameters
using LinearAlgebra: norm
using ..ParameterSets: P1
using ..References: k0, dD0, gamma_r, host_roots, green_1d_patch, kappa1
using ..FEBuilders: solve_configA, graded_mesh
using ..Observables: fit_waves, sample_plate, incoming_field

export L1_1_params, run_L1_1

"""
Case of L1-1: host (`T_ρ`, depth `H`), frequency `ω`, mesh
(`elements_per_wavelength`, `order`, vertical grading `max_ratio`) and layout in
plate wavelengths: the window spans ±`half_window` λ_p around the source, and
the probes start `probe_offset` λ_p from it. The buffer to the sponges is the
builder default (source ≥ 5/κ₁ from the sponges, `FEBuilders.configA_channel`):
with the sponges 11.5 m from the source at ω = 10 (P1, 1/κ₁ = 6.3 m), the
evanescent reflection left a mesh-independent error of 9e-6 in Re G.
"""
@with_kw struct L1_1_params
    T_ρ::Float64 = 0.0
    H::Float64 = 10.0
    ω::Float64 = 10.0
    elements_per_wavelength::Int = 20
    order::Int = 2
    max_ratio::Float64 = 1.2
    half_window::Float64 = 4.0
    probe_offset::Float64 = 1.0
    n_modes::Int = 6
end

"Run one L1-1 case; returns the raw FE measures and their references."
function run_L1_1(p::L1_1_params)
    @unpack T_ρ, H, ω, elements_per_wavelength, order, max_ratio, half_window, probe_offset, n_modes = p
    h = P1(T_ρ=T_ρ).host
    h = typeof(h)(h; H=H)
    F = 1.0
    kp = k0(h, ω)
    λ = 2π / kp
    mesh = graded_mesh(h, ω, elements_per_wavelength; order=order, max_ratio=max_ratio)
    sol = solve_configA(h, ω; window_length=2half_window * λ, source_x=half_window * λ, F=F, mesh=mesh)
    lay = sol.layout
    xs, b = lay.source_x, lay.source_width
    roots = host_roots(h, ω)

    # far field: discrete wavenumber and outgoing/incoming amplitudes
    right = fit_waves(sol.ηh, roots, xs + probe_offset * λ, lay.window[2]; x_ref=xs, n_modes=n_modes)
    left = fit_waves(sol.ηh, roots, lay.window[1], xs - probe_offset * λ; x_ref=xs, n_modes=n_modes)
    out_R, in_R = right.Aplus, right.Aminus
    out_L, in_L = left.Aminus, left.Aplus
    s = sin(kp * b / 2) / (kp * b / 2)
    γ_R = (out_R - in_L) / (im * s * F)
    γ_L = (out_L - in_R) / (im * s * F)

    # near field on the right of the patch, incoming sponge waves removed
    # (the source is transparent to them)
    kh = right.k
    xn = collect(range(xs + b / 2, xs + probe_offset * λ; length=401))
    incoming_prop(x) = in_L * exp(im * kh * (x - xs)) + in_R * exp(-im * kh * (x - xs))
    ηraw = sample_plate(sol.ηh, xn)
    ηn = ηraw .- incoming_prop.(xn)
    ηe = ηraw .- incoming_field.(Ref(left), Ref(right), xn)
    ηref = [F * green_1d_patch(h, max(x - xs, b / 2), ω, b; roots=roots) for x in xn]   # round-off at the edge
    rel(a, b) = norm(a - b) / norm(b)

    return Dict{String,Any}(
        "k0" => kp, "k_h_R" => right.k, "k_h_L" => left.k,
        "gamma_r" => gamma_r(h, ω; dim=2), "gamma_far_R" => γ_R, "gamma_far_L" => γ_L,
        "x_near" => xn .- xs, "eta_near_h" => ηn, "eta_near_ref" => ηref,
        "err_near" => rel(ηn, ηref), "err_near_re" => rel(real(ηn), real(ηref)),
        "err_near_im" => rel(imag(ηn), imag(ηref)),
        "err_near_re_evan" => rel(real(ηe), real(ηref)), "err_near_im_evan" => rel(imag(ηe), imag(ηref)),
        "eta_edge_h" => first(ηn), "eta_edge_ref" => first(ηref),
        "sponge_reflection" => max(abs(in_R / out_R), abs(in_L / out_L)),
        "fit_residual" => max(right.residual, left.residual),
        "lambda_p" => λ, "b" => b, "element_size" => lay.element_size, "nx" => lay.nx,
        "sponge_distance" => xs - lay.sponge_length, "kappa1" => kappa1(h, ω),
        "ny" => mesh.ny, "vertical_ratio" => mesh.vertical_ratio,
        "ndofs" => length(sol.result.solution.free_values),
    )
end

end # module L1_1_BareHostDispersion
