"""
Unit tests for the verification-suite workflow layer (src/common/workflow.jl).
"""

using Test
using HydroElasticFEM_xxx2026
const HX = HydroElasticFEM_xxx2026

Base.@kwdef struct DummyCase
    name::String = "dummy"
    n::Int = 4
    ω::Float64 = 1.5
end

@testset "Workflow: output paths" begin
    mktempdir() do tmp
        with_output_root(tmp) do
            @test output_root() == abspath(tmp)
            @test sims_path("L1-5") == joinpath(abspath(tmp), "data", "sims", "L1-5")
            @test exp_pro_path("L2-3", "a.csv") == joinpath(abspath(tmp), "data", "exp_pro", "L2-3", "a.csv")
            @test plots_path("L4-4", "f.pdf") == joinpath(abspath(tmp), "plots", "L4", "L4-4", "f.pdf")
        end
        @test output_root() == HX.REPO_ROOT
    end
end

@testset "Workflow: produce_or_load_case with provenance" begin
    mktempdir() do tmp
        with_output_root(tmp) do
            calls = Ref(0)
            run(case) = (calls[] += 1; Dict{String,Any}("result" => case.n * case.ω))
            data, file = produce_or_load_case(run, DummyCase(), "L0-1")
            @test calls[] == 1
            @test isfile(file)
            @test startswith(file, sims_path("L0-1"))
            @test data["result"] == 6.0
            for key in ("test_id", "matrix_version", "hefem_version", "hefem_commit",
                        "julia_version", "nthreads", "wall_time_s", "gitcommit")
                @test haskey(data, key)
            end
            @test data["test_id"] == "L0-1"

            data2, file2 = produce_or_load_case(run, DummyCase(), "L0-1")
            @test calls[] == 1          # loaded, not recomputed
            @test file2 == file
            @test data2["result"] == 6.0

            produce_or_load_case(run, DummyCase(), "L0-1"; force=true)
            @test calls[] == 2
            produce_or_load_case(run, DummyCase(n=8), "L0-1")
            @test calls[] == 3          # distinct case → distinct record
            @test length(readdir(sims_path("L0-1"))) == 2

            @test_throws ErrorException produce_or_load_case(c -> 1.0, DummyCase(n=9), "L0-1")
        end
    end
end

@testset "Workflow: acceptance checks" begin
    c = HX.AcceptanceCheck("L1-5", "omega_X"; value=9.6461, reference=9.646, tolerance=1e-4)
    @test c.error ≈ 1e-4 / 9.646
    @test c.passed
    c0 = HX.AcceptanceCheck("L2-2", "woodbury"; value=3e-13, reference=0.0, tolerance=1e-10, kind="I")
    @test c0.error == 3e-13 && c0.passed
    @test !HX.AcceptanceCheck("L1-1", "k"; value=1.1, reference=1.0, tolerance=1e-3).passed
    @test_throws ArgumentError HX.AcceptanceCheck("L1-1", "k"; value=1.0, reference=1.0, tolerance=1.0, kind="X")

    mktempdir() do tmp
        with_output_root(tmp) do
            checks = AcceptanceCheck[]
            add_check!(checks, "L1-5", "omega_X"; value=9.6461, reference=9.646, tolerance=1e-4)
            add_check!(checks, "L1-5", "omega_G"; value=9.99, reference=9.986, tolerance=1e-4)
            path = write_checks("L1-5", checks)
            lines = readlines(path)
            @test length(lines) == 3
            @test startswith(lines[1], "test_id,name,kind,value,reference,error,tolerance,passed")
            @test endswith(lines[2], "true") && endswith(lines[3], "false")
        end
    end
end

@testset "Workflow: figures and LaTeX tables" begin
    mktempdir() do tmp
        with_output_root(tmp) do
            df = HX.DataFrame(case=["T_rho=0", "T_rho=9.8"], omega_X=[9.646, 9.916], passed=[true, false])
            path = write_latex_table("L1-5", "band_edges", df; caption="Band edges", label="tab:L1-5",
                                     formatters=Dict("omega_X" => x -> HX.@sprintf("%.3f", x)))
            tex = read(path, String)
            @test occursin("\\label{tab:L1-5}", tex)
            @test occursin("T\\_rho=0 & 9.646 & yes \\\\", tex)
            @test occursin("\\bottomrule", tex)

            fig = HX.CairoMakie.Figure()
            HX.CairoMakie.lines!(HX.CairoMakie.Axis(fig[1, 1]), 1:3, [1.0, 4.0, 9.0])
            files = save_figure(fig, "L0-3", "energy")
            @test isfile(files.pdf) && isfile(files.png)
            @test dirname(files.pdf) == plots_path("L0-3")
        end
    end
end

@testset "Workflow: test registry" begin
    @test_throws ArgumentError HX.register_test("X1-1", () -> nothing)
    received = Ref{Any}(nothing)
    entry(; kwargs...) = (received[] = values(kwargs); :ok)
    HX.register_test("L4-98", entry)
    HX.register_test("L4-99", entry)
    HX.register_test("L3-98", entry)
    try
        @test findfirst(==("L3-98"), list_tests()) < findfirst(==("L4-98"), list_tests())
        @test run_tests("L4-98"; test_suit=true) == Dict("L4-98" => :ok)
        @test received[].test_suit == true
        @test received[].make_plots == true
        @test keys(run_tests("L4")) ⊇ Set(["L4-98", "L4-99"])
        @test_throws ErrorException run_tests("L9-1")
    finally
        foreach(id -> delete!(HX.TEST_REGISTRY, id), ("L4-98", "L4-99", "L3-98"))
    end
end
