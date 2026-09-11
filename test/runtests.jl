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

@testset "Density-graded dispersion" begin
    plate, resonator = liu_2025_parameters()
    grading = DensityGradingParameters(resonator.n0, 1.0)
    @test isapprox(n_density_graded(0.0, grading), resonator.n0)
    @test isapprox(n_density_graded(1.0, grading), resonator.n0 / exp(1))
    @test isapprox(Mᵣ(resonator, 0.5), 5.0)
    branches = density_graded_dispersion(plate, resonator,
        n_density_graded(0.5, grading), 1.0)
    @test 0 < branches.lower < branches.upper
    table = density_graded_table(plate, resonator, grading,
        collect(range(-pi, pi; length=5)), [0.0, 0.5, 1.0])
    @test size(table.lower) == (3, 5)
    @test first(table.densities) > last(table.densities)
    @test all(table.local_frequency .== resonator.natural_frequency)
end

@testset "Frequency-graded dispersion" begin
    plate, resonator = liu_2025_parameters()
    grading = FrequencyGradingParameters(10.0, 1.0, 1.0)
    @test isapprox(ωᵣ_linear(0.0, grading), 10.0)
    @test isapprox(ωᵣ_linear(0.9, grading), 1.0)
    @test isapprox(ωᵣ_exponential(log(10.0), grading), 1.0)
    branches = graded_dispersion(plate, resonator, ωᵣ_linear(0.5, grading), 1.0)
    @test 0 < branches.lower < branches.upper
    table = frequency_graded_table(plate, resonator, grading,
        collect(range(-pi, pi; length=5)), [0.0, 0.5, 1.0])
    @test size(table.lower) == (3, 5)
    @test size(table.upper) == (3, 5)
end

@testset "Reproduction orchestration" begin
    plate, resonator = liu_2025_parameters()
    normalized_wave_numbers = collect(range(-pi, pi; length=5))
    table = dispersion_table(plate, resonator, normalized_wave_numbers)
    @test length(table.bare) == 5
    @test length(table.lower) == 5
    @test first(table.k_over_sqrt_n0) == -pi
    @test first(bandgap_limits(plate, resonator, normalized_wave_numbers)) >= 0
end