"""
Post-processing kernels shared by the verification tests: sampling FE fields,
error norms against references, and convergence-rate estimation.
"""
module Observables

using Gridap
using LinearAlgebra

export sample_plate, relative_error, observed_rates, richardson_extrapolate
export fit_waves, symmetric_scattering, patch_average, secant_root, half_power_width
export self_convergence, incoming_field

"Plate deflection `ηh` sampled at `xs` on the free-surface line z = `z`."
sample_plate(ηh, xs::AbstractVector{<:Real}; z::Real=0.0) =
    ComplexF64[ηh(Point(x, z)) for x in xs]

"Discrete relative 2-norm ‖a − b‖/‖b‖."
relative_error(a::AbstractVector, b::AbstractVector) = norm(a - b) / norm(b)

"""
    observed_rates(hs, errs)

Observed convergence orders log(eᵢ/eᵢ₊₁)/log(hᵢ/hᵢ₊₁) between consecutive levels.
"""
function observed_rates(hs::AbstractVector{<:Real}, errs::AbstractVector{<:Real})
    length(hs) == length(errs) || throw(DimensionMismatch("one error per mesh size"))
    return [log(errs[i] / errs[i + 1]) / log(hs[i] / hs[i + 1]) for i in 1:(length(hs) - 1)]
end

"""
    richardson_extrapolate(h1, f1, h2, f2, p)

Richardson extrapolation of a quantity computed at mesh sizes h1 > h2 with
convergence order p.
"""
richardson_extrapolate(h1, f1, h2, f2, p) = f2 + (f2 - f1) / ((h1 / h2)^p - 1)

"""
    self_convergence(hs, fs)

Observed order and Richardson extrapolation of a quantity `fs` (real or
complex) on the three finest mesh sizes `hs`, from the differences between
levels only: p = log(|f₂ − f₁|/|f₃ − f₂|)/log(h₂/h₃). Unlike a rate of
|f − f_ref|, this does not use the reference, so the extrapolated value can be
compared with it without circularity. Returns `(rate, value)`.
"""
function self_convergence(hs::AbstractVector{<:Real}, fs::AbstractVector)
    length(hs) >= 3 || throw(ArgumentError("self-convergence needs three levels"))
    (h1, h2, h3), (f1, f2, f3) = hs[(end - 2):end], fs[(end - 2):end]
    p = log(abs(f2 - f1) / abs(f3 - f2)) / log(h2 / h3)
    return (rate=p, value=richardson_extrapolate(h2, f2, h3, f3, p))
end

# ─── Wave amplitudes from FE fields ─────────────────────────────────────────

"""
    fit_waves(xs, ηs, roots; x_ref=0.0, n_modes=6, k=nothing)
    fit_waves(ηh, roots, xa, xb; npts=400, kwargs...)

Least-squares fit of a source-free plate field sampled at `xs` (or of `ηh` on
[xa, xb]) to the homogeneous solutions of the host,

    η(x) = A⁺e^{ik(x−x_ref)} + A⁻e^{−ik(x−x_ref)}
           + Σₙ cₙ⁺e^{ikₙ(x−x_a)} + cₙ⁻e^{ikₙ(x_b−x)},

with the `n_modes` non-propagating roots kₙ = `roots[2:n_modes+1]`
(Im kₙ > 0) decaying away from either end of the window. Fitting the
evanescent and flexural modes replaces the geometric separation of the probe
window from sources and scatterers. Keep `n_modes` small (the flexural pair
and a few evanescent modes): deeper modes are nearly zero across the window,
differ most between the FE and the continuum, and make the fit ill-conditioned
(with 10 modes, R of L1-3 jumped by 2e-5 for a 1e-9 change in ω). With `k = nothing` the propagating
wavenumber is fitted too (residual scanned on a grid within ±5 % of
k₀ = `roots[1]`, then golden-section search around the best grid point; the
residual has side minima on short windows), so that the discrete dispersion error does not leak into
the phases of A±; otherwise `k` is used as given.
Returns `(Aplus, Aminus, k, residual)` with residual = ‖fit − η‖/‖η‖, and
the evanescent part (`ks`, coefficients `cplus` of e^{ikₙ(x−x_a)} and `cminus`
of e^{ikₙ(x_b−x)}, window `xa`, `xb`, `x_ref`) used by `incoming_field`.
"""
function fit_waves(xs::AbstractVector{<:Real}, ηs::AbstractVector, roots::AbstractVector;
                   x_ref::Real=0.0, n_modes::Int=6, k::Union{Nothing,Real}=nothing)
    xa, xb = extrema(xs)
    ks = roots[2:min(n_modes + 1, length(roots))]
    η = ComplexF64.(ηs)
    function solve(kk)
        cols = [exp.(im * kk .* (xs .- x_ref)), exp.(-im * kk .* (xs .- x_ref))]
        for kn in ks
            push!(cols, exp.(im * kn .* (xs .- xa)), exp.(im * kn .* (xb .- xs)))
        end
        M = reduce(hcat, cols)
        c = M \ η
        return c, norm(M * c - η) / norm(η)
    end
    if k === nothing
        k0 = real(roots[1])
        grid = range(0.95k0, 1.05k0; length=201)
        i = argmin([last(solve(kk)) for kk in grid])
        k = _golden_min(kk -> last(solve(kk)), grid[max(i - 1, 1)], grid[min(i + 1, end)]; tol=1e-13k0)
    end
    c, res = solve(k)
    return (Aplus=c[1], Aminus=c[2], k=k, residual=res, x_ref=x_ref, xa=xa, xb=xb,
            ks=ks, cplus=c[3:2:end], cminus=c[4:2:end])
end

"Minimiser of a unimodal `f` on [a, b] by golden-section search."
function _golden_min(f, a::Real, b::Real; tol::Real=1e-12)
    φ = (sqrt(5) - 1) / 2
    c, d = b - φ * (b - a), a + φ * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol
        if fc < fd
            b, d, fd = d, c, fc
            c = b - φ * (b - a)
            fc = f(c)
        else
            a, c, fc = c, d, fd
            d = a + φ * (b - a)
            fd = f(d)
        end
    end
    return (a + b) / 2
end

fit_waves(ηh, roots::AbstractVector, xa::Real, xb::Real; npts::Int=400, kwargs...) =
    (xs = range(xa, xb; length=npts); fit_waves(xs, sample_plate(ηh, xs), roots; kwargs...))

"""
    incoming_field(left, right, x)

Field entering a transparent region (e.g. a source patch, no scatterer)
from the domain ends, evaluated at `x` between the two fitted windows
`left` (x < region) and `right` (x > region) of `fit_waves`: the propagating
waves travelling towards the region (sponge reflections) plus the evanescent
modes decaying away from the far ends of the windows (reflections of the
evanescent near field at the sponges and walls). Subtracting it from the FE
field leaves the response to the source alone in an unbounded host.
"""
function incoming_field(left, right, x::Real)
    z = left.Aplus * exp(im * left.k * (x - left.x_ref)) +
        right.Aminus * exp(-im * right.k * (x - right.x_ref))
    z += sum(left.cplus .* exp.(im .* left.ks .* (x - left.xa)); init=zero(ComplexF64))
    z += sum(right.cminus .* exp.(im .* right.ks .* (right.xb - x)); init=zero(ComplexF64))
    return z
end

"""
    symmetric_scattering(up, down)

Reflection and transmission of a mirror-symmetric scatterer from the wave
amplitudes fitted upstream (`up`) and downstream (`down`), both referenced to
its centre of symmetry (`x_ref` of `fit_waves`):

    down.Aplus  = T·up.Aplus + R·down.Aminus
    up.Aminus   = R·up.Aplus + T·down.Aminus

Waves entering from the downstream side (sponge reflections) are removed
exactly. Returns `(R, T)` for incidence from upstream.
"""
function symmetric_scattering(up, down)
    M = [up.Aplus down.Aminus; down.Aminus up.Aplus]
    T, R = M \ [down.Aplus, up.Aminus]
    return (R=R, T=T)
end

"""
    patch_average(ηh, x_lo, x_hi, nel; npts=8)

Mean of the plate deflection over [x_lo, x_hi], made of `nel` elements, with
`npts` Gauss–Legendre points per element (exact for polynomial degree
≤ 2npts − 1).
"""
function patch_average(ηh, x_lo::Real, x_hi::Real, nel::Int; npts::Int=8)
    ξ, w = _gauss_legendre(npts)
    he = (x_hi - x_lo) / nel
    s = zero(ComplexF64)
    for e in 0:(nel - 1), (ξq, wq) in zip(ξ, w)
        s += wq / 2 * ηh(Point(x_lo + (e + (ξq + 1) / 2) * he, 0.0))
    end
    return s / nel
end

"Gauss–Legendre nodes and weights on [−1, 1] (Golub–Welsch)."
function _gauss_legendre(n::Int)
    β = [i / sqrt(4i^2 - 1) for i in 1:(n - 1)]
    F = eigen(SymTridiagonal(zeros(n), β))
    return F.values, 2 .* F.vectors[1, :] .^ 2
end

# ─── Scalar post-processing ─────────────────────────────────────────────────

"""
    secant_root(f, x0, x1; xtol=1e-12, maxiter=30)

Secant iteration for a root of the real scalar `f`; one evaluation per
iteration, for observables that need an FE solve per evaluation.
Returns `(root, evaluations, converged)`.
"""
function secant_root(f, x0::Real, x1::Real; xtol::Real=1e-12, maxiter::Int=30)
    f0, f1 = f(x0), f(x1)
    n = 2
    while abs(x1 - x0) > xtol * abs(x1) && n < maxiter && f1 != f0
        x0, x1 = x1, x1 - f1 * (x1 - x0) / (f1 - f0)
        f0, f1 = f1, f(x1)
        n += 1
    end
    return (root=x1, evaluations=n, converged=abs(x1 - x0) <= xtol * abs(x1))
end

"""
    half_power_width(xs, ys)

Full width of the highest peak of `ys` at half its maximum, by linear
interpolation; `NaN` when the peak is not bracketed within `xs`.
"""
function half_power_width(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    i = argmax(ys)
    half = ys[i] / 2
    l = findlast(j -> ys[j] < half, 1:i)
    r = findfirst(j -> ys[j] < half, i:length(ys))
    (l === nothing || r === nothing) && return NaN
    r += i - 1
    cross(a, b) = xs[a] + (half - ys[a]) * (xs[b] - xs[a]) / (ys[b] - ys[a])
    return cross(r - 1, r) - cross(l, l + 1)
end

end # module Observables
