"""
L1-3 — Single line resonator (`eq:lat_detuning`, `eq:lat_renorm_single`,
`eq:lat_absorbed_1d`).

One resonator on the P1 plate (configuration A), driven by a wave from a patch
source upstream. R and T are taken from the waves fitted on probe windows on
either side of the resonator (`Observables.symmetric_scattering`, phase
reference at the resonator), so that sponge reflections do not enter. The
propagating wavenumber of these fits is the discrete one, k_h = k₀(ω)(1 + δ),
with the relative dispersion error δ measured once per case on the bare plate
at ω_r (long one-way window); fitting k on the short probe windows themselves
is ill-posed when the field is a standing wave.

For a line resonator, R = iγ_r/(x − iΓ) holds for the FE with its own discrete
G_h, so 1/R = −(Γ + ix)/γ_r: the renormalised resonance x = 0 is the root of
Im(1/R), and the half-power points |x| = Γ of both |R|² and 𝒜 are the roots of
Im(1/R) ± Re(1/R). These smooth functions are solved on the FE by secant
iteration, and compared with the same roots of the exact single-resonator
solution (`References.single_resonator_1d`). A root that does not converge
within the swept band ω_r(1 ± span/Q_rad) is reported as NaN (unresolved
resonance, e.g. P2 at 8 elements per wavelength).
"""
module L1_3_SingleLineResonator

using Parameters
using Roots: find_zero
using ..ParameterSets: P1, ResonatorSpec, resonator_from_frequency
using ..References: k0, single_resonator_1d, renormalised_frequency, critical_damping, Q_rad, host_roots
using ..FEBuilders: solve_configA, ladder_mesh
using ..Observables: fit_waves, symmetric_scattering, secant_root, sample_plate

export L1_3_params, run_L1_3

"""
Case of L1-3. The dashpot is `c_factor` times the critical value c_crit
(γ_d = γ_r at the renormalised frequency of the reference). `snap` places the
resonator on a mesh node or at an element midpoint. The mesh is the ladder
level `level` at ω_r. With `n_omega > 0` the response is also swept over
ω_r(1 ± `span`/Q_rad). Lengths are in plate wavelengths at ω_r.
"""
@with_kw struct L1_3_params
    T_ρ::Float64 = 0.0
    m_ρr::Float64 = 0.01
    ω_r::Float64 = 10.0
    c_factor::Float64 = 0.0
    snap::Symbol = :node
    level::Int = 1
    order::Int = 2
    n_omega::Int = 0
    span::Float64 = 5.0
    probe::Float64 = 3.0
    offset::Float64 = 1.0
    n_modes::Int = 6
end

"Resonator of the case and the reference values that define it."
function L1_3_resonator(p::L1_3_params)
    h = P1(T_ρ=p.T_ρ).host
    r0 = resonator_from_frequency(p.m_ρr, p.ω_r)
    ω̃ = renormalised_frequency(h, r0)
    c_crit = critical_damping(h, r0, ω̃)
    return h, ResonatorSpec(m_ρ=r0.m_ρ, k_ρ=r0.k_ρ, c_ρ=p.c_factor * c_crit), ω̃, c_crit
end

"""
FE scattering of the case at frequency ω: R (phase reference at the
resonator), T, 𝒜 = 1 − |R|² − |T|², attachment and resonator displacements.
"""
function L1_3_scattering(p::L1_3_params, h, r::ResonatorSpec, ω::Real; δ::Real=0.0)
    λ, x1, window_length, mesh = _L1_3_layout(p, h)
    sol = solve_configA(h, ω; window_length=window_length, source_x=0.25λ,
                        resonators=[r], resonator_xs=[x1], snap_resonators=p.snap, mesh=mesh,
                        design_ω=p.ω_r)
    xr = only(sol.layout.resonator_xs)
    roots = host_roots(h, ω)
    kh = real(roots[1]) * (1 + δ)
    up = fit_waves(sol.ηh, roots, xr - (p.offset + p.probe) * λ, xr - p.offset * λ;
                   x_ref=xr, n_modes=p.n_modes, k=kh)
    down = fit_waves(sol.ηh, roots, xr + p.offset * λ, xr + (p.offset + p.probe) * λ;
                     x_ref=xr, n_modes=p.n_modes, k=kh)
    s = symmetric_scattering(up, down)
    η1 = only(sample_plate(sol.ηh, [xr]))
    return (R=s.R, T=s.T, A=1 - abs2(s.R) - abs2(s.T), q=only(sol.q), η1=η1, Ainc=up.Aplus,
            fit_residual=max(up.residual, down.residual), element_size=sol.layout.element_size,
            offset_in_element=mod(xr / sol.layout.element_size, 1.0))
end

"Window geometry (offsets from the window start) and mesh of a case."
function _L1_3_layout(p::L1_3_params, h)
    λ = 2π / k0(h, p.ω_r)
    x1 = 0.25λ + (0.5 + p.probe + p.offset) * λ
    return λ, x1, x1 + (p.offset + p.probe + 0.25) * λ, ladder_mesh(h, p.ω_r, p.level; order=p.order)
end

"""
Relative dispersion error δ = k_h/k₀ − 1 at ω_r on the mesh of the case, from
the bare plate driven by the same source (wavenumber fitted on the one-way
field from one wavelength past the source to the end of the window).
"""
function L1_3_dispersion_error(p::L1_3_params, h)
    λ, _, window_length, mesh = _L1_3_layout(p, h)
    sol = solve_configA(h, p.ω_r; window_length=window_length, source_x=0.25λ, mesh=mesh, design_ω=p.ω_r)
    roots = host_roots(h, p.ω_r)
    fit = fit_waves(sol.ηh, roots, sol.layout.source_x + λ, sol.layout.window[2]; n_modes=p.n_modes)
    return fit.k / real(roots[1]) - 1
end

"Run one L1-3 case; returns peak, half-power and (optionally) sweep data."
function run_L1_3(p::L1_3_params)
    h, r, ω̃_ref, c_crit = L1_3_resonator(p)
    Q = Q_rad(h, r, p.ω_r)
    δ = L1_3_dispersion_error(p, h)
    cache = Dict{Float64,Any}()
    fe(ω) = get!(() -> L1_3_scattering(p, h, r, ω; δ=δ), cache, ω)
    invR(ω) = 1 / fe(ω).R
    invR_ref(ω) = 1 / single_resonator_1d(h, r, ω).R

    # renormalised resonance x = 0 and half-power points x = ±Γ; x decreases
    # with ω, so x = Γ (root of Im − Re of 1/R) lies below ω̃ and x = −Γ above
    Δ = 4p.ω_r / Q
    lower(f) = ω -> imag(f(ω)) - real(f(ω))
    upper(f) = ω -> imag(f(ω)) + real(f(ω))
    band = p.ω_r .* (1 - p.span / Q, 1 + p.span / Q)
    function fe_root(f, x0)
        s = secant_root(f, x0, x0 * (1 + 1e-3); xtol=1e-8)   # fit noise ≈ 1e-7 in f
        return s.converged && band[1] <= s.root <= band[2] ? s.root : NaN
    end
    ω̃_h = fe_root(ω -> imag(invR(ω)), ω̃_ref)
    ref_minus = find_zero(lower(invR_ref), (ω̃_ref - Δ, ω̃_ref))
    ref_plus = find_zero(upper(invR_ref), (ω̃_ref, ω̃_ref + Δ))
    ω_minus = fe_root(lower(invR), ref_minus)
    ω_plus = fe_root(upper(invR), ref_plus)

    peak, ref = fe(isnan(ω̃_h) ? ω̃_ref : ω̃_h), single_resonator_1d(h, r, ω̃_ref)
    κ = r.k_ρ - im * ω̃_h * r.c_ρ                 # resonator equation: q/η = κ/(κ − ω²m)
    data = Dict{String,Any}(
        "omega_peak_h" => ω̃_h, "omega_peak_ref" => ω̃_ref,
        "R_peak_h" => peak.R, "T_peak_h" => peak.T, "A_peak_h" => peak.A,
        "R_peak_ref" => ref.R, "T_peak_ref" => ref.T, "A_peak_ref" => ref.absorbed,
        "q_over_eta_peak_h" => peak.q / peak.η1,
        "q_over_eta_peak_ref" => κ / (κ - ω̃_h^2 * r.m_ρ),
        "omega_half_h" => [ω_minus, ω_plus], "omega_half_ref" => [ref_minus, ref_plus],
        "Q_rad" => Q, "c_crit" => c_crit, "c_ρr" => r.c_ρ,
        "fit_residual_peak" => peak.fit_residual, "element_size" => peak.element_size,
        "offset_in_element" => peak.offset_in_element, "lambda_p" => 2π / k0(h, p.ω_r),
        "fe_solves_roots" => length(cache), "dispersion_error" => δ,
    )
    if p.n_omega > 0
        ωs = collect(range(p.ω_r * (1 - p.span / Q), p.ω_r * (1 + p.span / Q); length=p.n_omega))
        sweep = fe.(ωs)
        refs = single_resonator_1d.(Ref(h), Ref(r), ωs)
        merge!(data, Dict{String,Any}(
            "omega" => ωs,
            "R_h" => getfield.(sweep, :R), "T_h" => getfield.(sweep, :T), "A_h" => getfield.(sweep, :A),
            "R_ref" => getfield.(refs, :R), "T_ref" => getfield.(refs, :T), "A_ref" => getfield.(refs, :absorbed),
            "fit_residual_sweep" => maximum(getfield.(sweep, :fit_residual)),
        ))
    end
    return data
end

end # module L1_3_SingleLineResonator
