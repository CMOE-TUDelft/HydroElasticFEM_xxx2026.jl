"""
Smoke test of every registered verification test: the warm-up case
(`test_suit=true`) runs end to end, including CSV, checks, figures and tables,
in a temporary output root.
"""

using Test
using HydroElasticFEM_xxx2026

@testset "Smoke: $(id)" for id in list_tests()
    mktempdir() do tmp
        with_output_root(tmp) do
            result = run_tests(id; test_suit=true)[id]
            @test result.data isa HydroElasticFEM_xxx2026.DataFrame
            @test HydroElasticFEM_xxx2026.nrow(result.data) >= 1
            @test isfile(exp_pro_path(id, "checks.csv"))
            @test !isempty(readdir(plots_path(id)))
        end
    end
end
