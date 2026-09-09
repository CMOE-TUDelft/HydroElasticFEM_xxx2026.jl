using Test
using HydroElasticFEM_xxx2026

@testset "Appendix-A explicit dispersion" begin
    plate, resonator = liu_2025_parameters()
    k = 1.0
    bare = bare_dispersion(plate, k)
    branches = metaplate_dispersion(plate, resonator, k)

    @test isfinite(bare)
    @test 0 < branches.lower < branches.upper

    # With vanishing resonator density, the lower branch is the bare branch
    zero_resonator = ResonatorParameters(0.0, resonator.resonator_mass,
        resonator.natural_frequency)
    zero_branches = metaplate_dispersion(plate, zero_resonator, k)
    @test isapprox(zero_branches.lower, bare; rtol=1e-12)
    @test isapprox(zero_branches.upper, resonator.natural_frequency; rtol=1e-12)

    pretensioned, _ = liu_2025_parameters(tension=10.0)
    @test bare_dispersion(pretensioned, k) > bare
    @test_throws DomainError bare_dispersion(plate, 0.0)
end