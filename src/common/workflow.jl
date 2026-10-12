"""
Workflow utilities shared by every verification test: output locations,
DrWatson bookkeeping with provenance, acceptance checks, figure/table export
and the `run_tests` registry.

Output layout (relative to `output_root()`, the repository root by default):

    data/sims/<test-id>/        raw `produce_or_load` records (JLD2, git-ignored)
    data/exp_pro/<test-id>/     processed CSV, including `checks.csv` (committed)
    plots/<level>/<test-id>/    paper figures (.pdf, .png) and tables (.tex)
"""

import HydroElasticFEM
import TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"Version of the verification test matrix the cases are transcribed from."
const MATRIX_VERSION = "2026-10-12"

const _OUTPUT_ROOT = Ref{String}(REPO_ROOT)

# ─── Output locations ───────────────────────────────────────────────────────

"Root directory under which `data/` and `plots/` are written."
output_root() = _OUTPUT_ROOT[]

"""
    with_output_root(f, dir)

Run `f()` with all outputs redirected to `dir` (used by the smoke tests).
"""
function with_output_root(f, dir::AbstractString)
    previous = _OUTPUT_ROOT[]
    _OUTPUT_ROOT[] = abspath(dir)
    try
        return f()
    finally
        _OUTPUT_ROOT[] = previous
    end
end

"Level of a test id, e.g. `\"L1\"` for `\"L1-5\"`."
test_level(test_id::AbstractString) = first(split(test_id, '-'))

sims_path(test_id::AbstractString, args...) =
    joinpath(output_root(), "data", "sims", test_id, args...)
exp_pro_path(test_id::AbstractString, args...) =
    joinpath(output_root(), "data", "exp_pro", test_id, args...)
plots_path(test_id::AbstractString, args...) =
    joinpath(output_root(), "plots", test_level(test_id), test_id, args...)

# ─── Provenance and DrWatson bookkeeping ────────────────────────────────────

function _git_commit(path::AbstractString)
    try
        return readchomp(pipeline(`git -C $path rev-parse HEAD`; stderr=devnull))
    catch
        return "unknown"
    end
end

"""
    _hefem_pin()

HydroElasticFEM source as pinned in this repository's Manifest: repository
URL or dev path, tracked revision and git tree hash.
"""
function _hefem_pin()
    manifest = joinpath(REPO_ROOT, "Manifest.toml")
    isfile(manifest) || return (source="unknown", rev="unknown", tree="unknown")
    entry = only(TOML.parsefile(manifest)["deps"]["HydroElasticFEM"])
    return (source=get(entry, "repo-url", get(entry, "path", "registry")),
            rev=get(entry, "repo-rev", "unknown"),
            tree=get(entry, "git-tree-sha1", "unknown"))
end

"""
    provenance(test_id)

Metadata stored with every simulation record, in addition to the git commit
and patch of this repository added by DrWatson (`tag=true`).
`hefem_commit` is only known when HydroElasticFEM is a `dev` checkout;
otherwise the Manifest pin (`hefem_source`, `hefem_rev`, `hefem_tree`)
identifies the code.
"""
function provenance(test_id::AbstractString)
    pin = _hefem_pin()
    return Dict{String,Any}(
        "test_id"        => String(test_id),
        "matrix_version" => MATRIX_VERSION,
        "hefem_version"  => string(pkgversion(HydroElasticFEM)),
        "hefem_source"   => pin.source,
        "hefem_rev"      => pin.rev,
        "hefem_tree"     => pin.tree,
        "hefem_commit"   => _git_commit(pkgdir(HydroElasticFEM)),
        "julia_version"  => string(VERSION),
        "nthreads"       => Threads.nthreads(),
    )
end

"""
    produce_or_load_case(f, case, test_id; force=false, verbose=false, kwargs...)

Load the record of `case` from `sims_path(test_id)` or compute it with
`f(case)::Dict{String,Any}` and save it. Provenance, wall time and allocated
bytes are added to the record. Remaining `kwargs` are forwarded to
`DrWatson.savename` (e.g. `digits`, `accesses`).

Parameter structs whose array-valued fields matter must also carry a scalar
field (e.g. `layout_hash::String`) so that distinct cases get distinct names.
"""
function produce_or_load_case(f, case, test_id::AbstractString;
                              force::Bool=false, verbose::Bool=false, kwargs...)
    path = sims_path(test_id)
    mkpath(path)
    function compute(c)
        stats = @timed f(c)
        data = stats.value
        data isa AbstractDict || error("Test $(test_id): run function must return a Dict, got $(typeof(data))")
        merge!(data, provenance(test_id))
        data["wall_time_s"] = stats.time
        data["allocated_bytes"] = stats.bytes
        return data
    end
    return produce_or_load(compute, case, path; force=force, verbose=verbose,
                           tag=true, gitpath=REPO_ROOT, storepatch=true, kwargs...)
end

# ─── Acceptance checks ──────────────────────────────────────────────────────

"""
One pass/fail criterion of the verification test matrix.

`kind` is `"I"` (identity, round-off), `"D"` (discretisation) or
`"A"` (asymptotic). `error` is relative to `reference` unless the reference
is zero, in which case it is absolute.
"""
struct AcceptanceCheck
    test_id::String
    name::String
    kind::String
    value::Float64
    reference::Float64
    error::Float64
    tolerance::Float64
    passed::Bool
end

function AcceptanceCheck(test_id, name; value, reference, tolerance, kind="D")
    kind in ("I", "D", "A") || throw(ArgumentError("kind must be I, D or A, got $(kind)"))
    err = iszero(reference) ? abs(value) : abs(value - reference) / abs(reference)
    return AcceptanceCheck(test_id, name, kind, value, reference, err, tolerance, err <= tolerance)
end

"Append an `AcceptanceCheck` to `checks` and return it."
function add_check!(checks::Vector{AcceptanceCheck}, test_id, name; kwargs...)
    check = AcceptanceCheck(test_id, name; kwargs...)
    push!(checks, check)
    return check
end

"Write the checks of `test_id` to `data/exp_pro/<test-id>/checks.csv`."
function write_checks(test_id::AbstractString, checks::Vector{AcceptanceCheck})
    path = exp_pro_path(test_id, "checks.csv")
    mkpath(dirname(path))
    fields = fieldnames(AcceptanceCheck)
    rows = [NamedTuple{fields}(getfield.(Ref(c), fields)) for c in checks]
    CSV.write(path, DataFrame(rows))
    return path
end

# ─── Figures and tables ─────────────────────────────────────────────────────

"Save `fig` as `<name>.pdf` and `<name>.png` in `plots_path(test_id)`."
function save_figure(fig, test_id::AbstractString, name::AbstractString)
    dir = plots_path(test_id)
    mkpath(dir)
    pdf = joinpath(dir, name * ".pdf")
    png = joinpath(dir, name * ".png")
    save(pdf, fig)
    save(png, fig, px_per_unit=300 / 96)
    return (pdf=pdf, png=png)
end

_latex_cell(x::AbstractString) = replace(x, "_" => "\\_", "%" => "\\%", "&" => "\\&")
_latex_cell(x::Integer) = string(x)
_latex_cell(x::Bool) = x ? "yes" : "no"
_latex_cell(x::Real) = @sprintf("%.6g", x)
_latex_cell(x) = _latex_cell(string(x))

"""
    write_latex_table(test_id, name, df; caption="", label="", formatters=Dict())

Write `df` as a booktabs `table` environment to `plots_path(test_id, name * ".tex")`.
`formatters` maps column names to functions returning the cell string.
"""
function write_latex_table(test_id::AbstractString, name::AbstractString, df::DataFrame;
                           caption::AbstractString="", label::AbstractString="",
                           formatters::AbstractDict=Dict{String,Function}())
    path = plots_path(test_id, name * ".tex")
    mkpath(dirname(path))
    cols = names(df)
    open(path, "w") do io
        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        isempty(caption) || println(io, "\\caption{", caption, "}")
        isempty(label) || println(io, "\\label{", label, "}")
        println(io, "\\begin{tabular}{", repeat("l", length(cols)), "}")
        println(io, "\\toprule")
        println(io, join(_latex_cell.(cols), " & "), " \\\\")
        println(io, "\\midrule")
        for row in eachrow(df)
            cells = [haskey(formatters, c) ? formatters[c](row[c]) : _latex_cell(row[c]) for c in cols]
            println(io, join(cells, " & "), " \\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")
    end
    return path
end

# ─── Test registry ──────────────────────────────────────────────────────────

const TEST_REGISTRY = Dict{String,Function}()
const _TEST_ID = r"^(L[0-4]|S3)-([0-9]+)$"

"""
    register_test(test_id, entry_point)

Called by each `scripts/*/<id>-<name>.jl` file to make its entry point
`entry_point(; force, make_plots, save_csv, verbose, vtk_output, test_suit)`
available through `run_tests`.
"""
function register_test(test_id::AbstractString, entry_point::Function)
    occursin(_TEST_ID, test_id) || throw(ArgumentError("Invalid test id $(test_id); expected e.g. L1-5 or S3-2"))
    TEST_REGISTRY[String(test_id)] = entry_point
    return nothing
end

function _test_sort_key(test_id::AbstractString)
    m = match(_TEST_ID, test_id)
    return (m.captures[1], parse(Int, m.captures[2]))
end

"Registered test ids in execution order."
list_tests() = sort!(collect(keys(TEST_REGISTRY)); by=_test_sort_key)

"""
    run_tests(name="all"; force=false, make_plots=true, save_csv=true,
              verbose=false, vtk_output=false, test_suit=false)

Run one test (`"L1-5"`), one level (`"L1"`) or every registered test (`"all"`).
With `test_suit=true` each test only runs its warm-up case (smoke testing).
Returns a `Dict` from test id to the value returned by its entry point.
"""
function run_tests(name::AbstractString="all";
                   force::Bool=false, make_plots::Bool=true, save_csv::Bool=true,
                   verbose::Bool=false, vtk_output::Bool=false, test_suit::Bool=false)
    ids = if name == "all"
        list_tests()
    elseif haskey(TEST_REGISTRY, name)
        [String(name)]
    else
        filter(id -> test_level(id) == name, list_tests())
    end
    isempty(ids) && name != "all" &&
        error("Unknown test or level: $(name). Registered tests: $(join(list_tests(), ", "))")
    results = Dict{String,Any}()
    for id in ids
        verbose && println("[run_tests] ", id)
        results[id] = TEST_REGISTRY[id](; force=force, make_plots=make_plots, save_csv=save_csv,
                                         verbose=verbose, vtk_output=vtk_output, test_suit=test_suit)
    end
    return results
end
