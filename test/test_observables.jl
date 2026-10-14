"""
Unit tests of the post-processing kernels (src/common/observables.jl) on
synthetic fields with known amplitudes.
"""

using Test
using LinearAlgebra
using HydroElasticFEM_xxx2026
using HydroElasticFEM_xxx2026.ParameterSets
using HydroElasticFEM_xxx2026.References
using HydroElasticFEM_xxx2026.Observables

@testset "Observables: modal wave fit" begin
    h, ω = P1().host, 10.0
    roots = host_roots(h, ω; n_evanescent=20)
    kp = real(roots[1])
    xa, xb, xr = 3.0, 12.0, 1.5
    Ap, Am = 0.7 - 0.2im, 0.1 + 0.3im
    # waves plus the slowest evanescent and the flexural modes from both ends
    field(x; k=kp) = Ap * exp(im * k * (x - xr)) + Am * exp(-im * k * (x - xr)) +
                     0.02exp(im * roots[4] * (x - xa)) + 0.01im * exp(im * roots[4] * (xb - x)) +
                     0.05exp(im * roots[2] * (x - xa)) + 0.03exp(im * roots[3] * (xb - x))
    xs = range(xa, xb; length=300)
    fit = fit_waves(xs, field.(xs), roots; x_ref=xr, n_modes=6, k=kp)
    @test fit.Aplus ≈ Ap rtol = 1e-10
    @test fit.Aminus ≈ Am rtol = 1e-10
    @test fit.residual < 1e-10
    # the propagating wavenumber is recovered when it is fitted
    kh = 0.997kp
    fitk = fit_waves(xs, field.(xs; k=kh), roots; x_ref=xr, n_modes=6)
    @test fitk.k ≈ kh rtol = 1e-9
    @test fitk.Aplus ≈ Ap rtol = 1e-7
    # without the evanescent modes the fit is visibly worse
    @test fit_waves(xs, field.(xs), roots; x_ref=xr, n_modes=0, k=kp).residual > 1e-3
end

@testset "Observables: symmetric scattering de-embedding" begin
    R, T = 0.3 + 0.4im, -0.5 + 0.2im
    a_in, b_in = 1.0 + 0.0im, 0.15 - 0.05im      # incident from upstream and downstream
    up = (Aplus=a_in, Aminus=R * a_in + T * b_in)
    down = (Aplus=T * a_in + R * b_in, Aminus=b_in)
    s = symmetric_scattering(up, down)
    @test s.R ≈ R
    @test s.T ≈ T
end

@testset "Observables: scalar helpers" begin
    @test secant_root(x -> x^3 - 2, 1.0, 1.5).root ≈ cbrt(2) rtol = 1e-12
    ωs = range(9.0, 11.0; length=4001)
    lorentz = @. 1 / (1 + ((ωs - 10.0) / 0.1)^2)
    @test half_power_width(ωs, lorentz) ≈ 0.2 rtol = 1e-4
    @test isnan(half_power_width(ωs, collect(ωs)))
    # self-convergence recovers order and limit without the reference
    hs = [0.4, 0.2, 0.1]
    sc = self_convergence(hs, @. 3.0 + 0.7hs^3 + 0.01hs^5)
    @test sc.rate ≈ 3 atol = 0.02
    @test sc.value ≈ 3.0 atol = 1e-5
    @test self_convergence(hs, @. (1 + 2im) + (0.5 - im) * hs^2).value ≈ 1 + 2im
    ξ, w = Observables._gauss_legendre(5)
    @test sum(w) ≈ 2
    @test sum(w .* ξ .^ 8) ≈ 2 / 9
end

@testset "Workflow: rate checks" begin
    checks = AcceptanceCheck[]
    @test add_rate_check!(checks, "L1-2", "fast"; rate=4.5, theoretical=4).passed
    @test add_rate_check!(checks, "L1-2", "slack"; rate=3.75, theoretical=4).passed
    @test !add_rate_check!(checks, "L1-2", "slow"; rate=3.5, theoretical=4).passed
end
