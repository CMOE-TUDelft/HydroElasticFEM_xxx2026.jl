"""
End-to-end tests of the configuration-A FE builder against the exact
references: patch-source Green function (bare plate) and Foldy solution
(line resonators), P1 host at ω ≈ ω_r.
"""

using Test
using HydroElasticFEM_xxx2026.ParameterSets
using HydroElasticFEM_xxx2026.References
using HydroElasticFEM_xxx2026.FEBuilders
using HydroElasticFEM_xxx2026.Observables

"FE plate field and Foldy reference (with attachment strength `αfun`) on the window."
function configA_vs_reference(h, ω, rs, offs; αfun=alpha, mesh=MeshSpec(elements_per_wavelength=24, ny=16))
    roots = host_roots(h, ω)
    sol = solve_configA(h, ω; window_length=8.0, source_x=0.5, resonators=rs, resonator_xs=offs, mesh=mesh)
    lay = sol.layout
    G(x) = green_1d(h, x, ω; roots=roots)
    Gb(x) = green_1d_patch(h, x, ω, lay.source_width; roots=roots)
    N = length(rs)
    α = [αfun(r, ω) for r in rs]
    Gm = [G(lay.resonator_xs[i] - lay.resonator_xs[j]) for i in 1:N, j in 1:N]
    ηatt = N == 0 ? ComplexF64[] :
           (I - Gm * Diagonal(α)) \ [Gb(x - lay.source_x) for x in lay.resonator_xs]
    field(x) = Gb(x - lay.source_x) + sum((α[j] * ηatt[j] * G(x - lay.resonator_xs[j]) for j in 1:N); init=0im)
    xs = range(lay.source_x + lay.source_width, lay.window[2]; length=81)
    return (field_error=relative_error(sample_plate(sol.ηh, xs), field.(xs)),
            attachment_error=N == 0 ? 0.0 : relative_error(sample_plate(sol.ηh, lay.resonator_xs), ηatt),
            sol=sol, ηatt=ηatt)
end

"Attachment strength with the plate-side dashpot force reversed (HydroElasticFEM before PR #75)."
alpha_reversed_dashpot(r::ResonatorSpec, ω) =
    ω^2 * r.m_ρ * (r.k_ρ + im * ω * r.c_ρ) / (r.k_ρ - im * ω * r.c_ρ - ω^2 * r.m_ρ)

using LinearAlgebra

@testset "Configuration A: bare plate vs patch Green function" begin
    h = P1().host
    c = configA_vs_reference(h, 10.0, ResonatorSpec[], Float64[])
    @test c.field_error < 2e-2                  # P2, 24 el/λ_p: ≈ 1 % (study in fe_builders.jl)
    lay = c.sol.layout
    # buffer: the source is ≥ 5/κ₁ from the sponges (and the buffer ≥ λ_p)
    @test lay.window[1] - lay.sponge_length >= lay.λ_p - 1e-12
    @test lay.source_x - lay.sponge_length >= 5 / kappa1(h, 10.0) - lay.element_size
    @test isapprox(lay.source_x / lay.element_size, round(lay.source_x / lay.element_size); atol=1e-9)
end

@testset "Configuration A: undamped resonators vs Foldy" begin
    ps = P1()
    h, r = ps.host, ps.resonator
    rs = [resonator_from_frequency(0.01, 10.0 + 0.2sin(3i)) for i in 1:3]
    c = configA_vs_reference(h, 9.8, rs, [2.0, 3.1, 4.4])
    @test c.field_error < 2e-2
    @test c.attachment_error < 2e-2
    # resonator displacement from its own equation of motion
    q_ref = [r.k_ρ / (r.k_ρ - 9.8^2 * r.m_ρ) for r in rs] .* c.ηatt
    @test relative_error(c.sol.q, q_ref) < 2e-2
end

@testset "Configuration A: damped resonator vs Foldy" begin
    h, r = P1().host, with_damping_ratio(P1().resonator, 0.05)
    mesh = MeshSpec(elements_per_wavelength=24, ny=16, order=3)
    paper = configA_vs_reference(h, 9.8, [r], [3.0]; mesh=mesh)
    reversed = configA_vs_reference(h, 9.8, [r], [3.0]; αfun=alpha_reversed_dashpot, mesh=mesh)
    @test paper.field_error < 5e-3
    @test paper.attachment_error < 5e-3
    # guards against the plate-side dashpot sign error fixed in HydroElasticFEM#75
    @test reversed.field_error > 5paper.field_error
end

@testset "Configuration A: mesh ladder and resonator snapping" begin
    h = P1().host
    m = ladder_mesh(h, 10.0, 1)
    @test m.elements_per_wavelength == 16 && m.ny == 16
    r = m.vertical_ratio
    @test h.H * (r - 1) / (r^m.ny - 1) ≈ 2π / k0(h, 10.0) / 32      # top cell = half the element size
    @test ladder_mesh(P1().host, 10.0, 0; order=3).order == 3
    frac(x) = abs(x - round(x))
    for (snap, expected) in ((:node, 0.0), (:midelement, 0.5))
        _, lay = configA_channel(h, 10.0; window_length=5.0, source_x=0.5,
                                 resonators=[P1().resonator], resonator_xs=[2.123], snap_resonators=snap)
        @test frac(only(lay.resonator_xs) / lay.element_size) ≈ expected atol = 1e-9
    end
    # the domain is sized at design_ω, independently of ω
    _, l1 = configA_channel(h, 9.5; window_length=5.0, source_x=0.5, design_ω=10.0)
    _, l2 = configA_channel(h, 10.5; window_length=5.0, source_x=0.5, design_ω=10.0)
    @test l1.L == l2.L && l1.nx == l2.nx
end
