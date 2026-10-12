"""
Post-processing kernels shared by the verification tests: sampling FE fields,
error norms against references, and convergence-rate estimation.
"""
module Observables

using Gridap
using LinearAlgebra

export sample_plate, relative_error, observed_rates, richardson_extrapolate

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

end # module Observables
