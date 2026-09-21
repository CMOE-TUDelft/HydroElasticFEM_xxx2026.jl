# Frequency-domain response of a floating metaplate with a DENSITY-graded
# resonator POPULATION, following the theory notes (Locally Resonant
# Hydroelasticity project), subsubsection "Density-graded resonator
# population":
#
#     ωᵣ(x) = ωᵣ  (fixed)          n(x) = n₀ e^{-x/L},  x ∈ [0,L]
#     M_r(x) = n(x) m_{ρ,r}
#
# i.e. the INDIVIDUAL resonator's own frequency never changes; what varies
# along the array is how many resonators per unit length are present. This
# is the density-grading counterpart of the attached frequency-graded
# script, which instead grades Kᵣ(i) (hence ωᵣ(i)) at UNIFORM spacing.
#
# DISCRETE REALISATION AND MESH-NODE ENFORCEMENT
# -----------------------------------------------
# A continuum number density n(x) has no unique discrete realisation, but
# the natural (and safest) one, re-using this codebase's own conventions,
# is:
#   1. Keep the SAME candidate lattice as the uniform/frequency-graded
#      arrays: Nᵣ_MAX candidate unit cells of width Lᵣ, cell centres at
#           ζ_i = Lᵣ/2 + (i-1)Lᵣ ,   i = 1,...,Nᵣ_MAX,
#      which the ORIGINAL mesh-sizing block (h = min(λ/20, Lᵣ/2), then
#      nxs chosen so Lᵣ/2 divides an integer number of h-elements) already
#      guarantees sit exactly on mesh nodes, x₀-relative.
#   2. "Populate" only a SUBSET of these candidate cells, chosen by a
#      deterministic error-diffusion / Bresenham-style thinning of the
#      target occupation probability p_i = n(ζ_i)/n₀ = e^{-ζ_i/L}. This is
#      the standard way to turn a spatially varying density into a point
#      pattern without introducing new candidate positions.
# Because every populated resonator is literally one of the untouched
# candidate cell-centres, it is on a mesh node BY CONSTRUCTION — no
# non-uniform re-meshing, no snapping/rounding, and no collision risk
# (populated cells are a subset of distinct integers, so positions are
# automatically distinct). We additionally assert this alignment explicitly
# below so the guarantee is checked, not just assumed, if you ever change
# Nᵣ_MAX, Lᵣ or the resolution rule.
#
# Mᵣ and Kᵣ are IDENTICAL for every populated resonator (both fixed, so
# ωᵣ = sqrt(Kᵣ/Mᵣ) is the same everywhere) — only the population/spacing
# encodes the grading, per the theory's "hold ωᵣ(x)=ωᵣ fixed" statement.
#
# ASSUMPTIONS YOU SHOULD CHECK / RETUNE:
#   - ω_r = 9.5 rad/s is a placeholder (midpoint of the frequency-graded
#     script's [8,11] rad/s band). Replace with your actual resonator
#     design frequency.
#   - L_grade = Ls, i.e. the grading (decay) length spans the full plate,
#     matching the theory figure caption (density falls from n/n₀=1 to
#     ≈0.37 = e^{-1} over x/L = 0 → 1).
#   - Everything else (tank geometry, plate/fluid properties, probes,
#     vertical mesh grading) is copied unchanged from the attached
#     uniform/frequency-graded scripts for direct comparability.
#
# NOTE: this script has not been executed — HydroElasticFEM and Julia are
# not available in the environment that produced it. It is written to
# mirror the attached, working script as closely as possible; please
# sanity-check the resonator_array call signature against your installed
# package version (see the comment at that call site).
#
# PARALLELISED WITH Distributed.jl, exactly as in the attached script: no
# shared memory across workers, so everything main() needs is pushed via
# @everywhere and pmap looks it up as a plain top-level (Main) name.

using Distributed

const N_WORKERS = max(Sys.CPU_THREADS - 1, 1)
if nprocs() == 1
  addprocs(N_WORKERS)
end

@everywhere begin
  using HydroElasticFEM
  using Gridap
  using LinearAlgebra
  import WaveSpec as WS
  import HydroElasticFEM.Geometry as G
  import HydroElasticFEM.Physics as P
  import HydroElasticFEM.Simulation as S
  import HydroElasticFEM.ParameterHandler as PH
  import HydroElasticFEM.PostProcessing as PP
  using HydroElasticFEM: map_vertical_GP_for_const_dep

  BLAS.set_num_threads(1)

  Ls = 30.0
  H = 10.0
  g = WS.PhysicalConstants.g
  mᵨ = 0.01
  EIᵨ = 0.05*g
  Mᵣ = 10.0            # resonator mass, CONSTANT across the population
  ρw = 1000.0
  Lᵣ = 1.0             # candidate unit-cell width (same lattice as the
                        # uniform / frequency-graded arrays)
  η0 = 0.01

  # ---------------------------------------------------------------------
  # Fixed resonator frequency (§ density-graded resonator population).
  # See the "ASSUMPTIONS" note at the top of the file.
  # ---------------------------------------------------------------------
  const ω_r = 9.5                 # rad/s, SAME for every populated resonator
  const Kᵣ  = Mᵣ * ω_r^2          # constant stiffness matching ω_r

  # ---------------------------------------------------------------------
  # Exponential density grading, n(x) = n₀ e^{-x/L_grade}, over the full
  # candidate lattice of Nᵣ_MAX = Ls/Lᵣ cells (i.e. n₀ = 1/Lᵣ is the fully-
  # populated reference density, one resonator per candidate cell).
  # ---------------------------------------------------------------------
  const Nᵣ_MAX  = Int(round(Ls / Lᵣ))     # = 30, matches the uniform array
  const L_grade = Ls                       # decay length spans the full plate

  # Deterministic error-diffusion thinning of the candidate lattice: keeps
  # cell i active whenever the running sum of its target occupation
  # probability p_i = n(ζ_i)/n₀ = e^{-ζ_i/L_grade} has accumulated to 1
  # "resonator's worth", then resets the accumulator. This is the discrete
  # analogue of ∫n(x)dx, so the realised population tracks the target
  # density profile without random seeding (fully reproducible) and never
  # revisits or duplicates a candidate index.
  function graded_occupation(nᵣ_max::Int, Lᵣ::Real, L_grade::Real)
    active = Int[]
    acc = 0.0
    for i in 1:nᵣ_max
      ζ_i = Lᵣ/2 + (i-1)*Lᵣ            # candidate cell-centre, x₀-relative
      p_i = exp(-ζ_i / L_grade)         # target n(ζ_i)/n₀
      acc += p_i
      if acc >= 1.0
        push!(active, i)
        acc -= 1.0
      end
    end
    return active
  end

  const ACTIVE_IDXS  = graded_occupation(Nᵣ_MAX, Lᵣ, L_grade)
  const Nᵣ_ACTIVE    = length(ACTIVE_IDXS)

  # Probe layout for R/T post-processing (unchanged from the attached script)
  const N_PROBES     = 4
  const PROBE_MARGIN = 1.5

  # Number of points sampled along the plate for the elevation-vs-frequency
  # map (see the field-sampling block in main()).
  const NX_FIELD_SAMPLE = 200

  function incident_wave(ω::Real, η0::Real)
    k = WS.AiryWaves.solve_wavenumber(ω, H)
    T = 2π / ω
    spec = WS.ContinuousSpectrums.RegularWave(2*η0, T)
    ds = WS.SpectralSpreading.DiscreteSpectralSpreading(spec; mess=false)
    spread = WS.AngularSpreading.DiscreteAngularSpreading(0.0)
    sea_state = WS.AiryWaves.AiryState(ds, spread, 1, 1, [ω], [k], [0.0], H, 1)
    ηin(x) = η0*exp(im*k*x[1])
    ϕin(x) = -im*(η0*ω/k)*(cosh(k*(x[2]+H)) / sinh(k*H))*exp(im*k*x[1])
    vin(x) = VectorValue((η0*ω)*(cosh(k*(x[2]+H)) / sinh(k*H))*exp(im*k*x[1]),
                         -im*(η0*ω)*(sinh(k*(x[2]+H)) / sinh(k*H))*exp(im*k*x[1]))
    return (; sea_state, ηin, ϕin, vin, k)
  end

  function domain_map(ny)
    # Vertical-only grading, identical to the uniform/frequency-graded
    # scripts. Density grading changes the resonator POPULATION, not the
    # plate/fluid mesh geometry, so x[1] is left untouched here.
    x -> VectorValue( x[1],
      map_vertical_GP_for_const_dep(x[2] - H, 1.08, ny, H; dbgmsg=false),
    )
  end

  function main(n::Int, ωinc::Real, vtk_output::Bool=false)

    inc_wave = incident_wave(ωinc, η0)
    f_in(x) = (inc_wave.vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * inc_wave.k * inc_wave.ϕin(x)

    λ = 2π / inc_wave.k
    println("Density-graded array — ωinc = $ωinc rad/s, k = $(inc_wave.k) rad/m, λ = $λ m")

    H_margin = min(H, λ)
    standoff = ceil(PROBE_MARGIN * H_margin + maximum(PP.suggest_probe_offsets(λ; n=N_PROBES)))
    # Mesh sizing identical to the uniform/frequency-graded scripts: h is
    # fine enough to resolve BOTH the wave (λ/20) and the finest possible
    # gap between active resonators, which is bounded below by Lᵣ (two
    # neighbouring candidate cells both active — happens near x=0 where
    # the occupation probability is close to 1). Reusing Lᵣ/2 here is
    # therefore still the right (and conservative) resolution requirement
    # even though the realised spacing elsewhere is larger due to thinning.
    h  = min(λ / 20, Lᵣ / 2)
    nxs = ceil(Int, (Lᵣ/2) / h)*ceil(Int, Ls/(Lᵣ/2))
    hs = Ls / nxs
    standoff_ = ceil(Int, standoff/hs)*hs
    L = 2*standoff_ + Ls
    x₀ = standoff_
    nx = Int(L / hs)
    println("Tank length L = $L m, nx = $nx, ny = $n, h = $h m, hs = $hs m")
    println("Populated $(Nᵣ_ACTIVE) of $(Nᵣ_MAX) candidate resonator cells " *
            "(density grading L=$(L_grade) m).")

    # --- Resonator positions: SAME candidate cell-centre formula as the
    # uniform/frequency-graded scripts, evaluated only at the populated
    # (ACTIVE_IDXS) indices. -----------------------------------------
    resonator_x = [x₀ + Lᵣ/2 + (i-1)*Lᵣ for i in ACTIVE_IDXS]

    # Explicit mesh-node-alignment check (enforcement, not just
    # assumption): every resonator offset from x₀ must be an integer
    # multiple of hs, since hs was constructed above to divide Lᵣ/2
    # exactly.
    for (k, i) in enumerate(ACTIVE_IDXS)
      offset = Lᵣ/2 + (i-1)*Lᵣ
      n_hs = offset / hs
      @assert abs(n_hs - round(n_hs)) < 1e-6 "Resonator $k (candidate cell $i, offset=$offset m) is not aligned to a mesh node (hs=$hs m) — check that Lᵣ/2 still divides hs exactly."
    end

    tank = G.TankDomain(
        L = L, H = H, nx = nx, ny = n,
        map=domain_map(n),
        is_periodic = (false, false),
        structure_domains = [G.StructureDomain(L=Ls, x₀=[x₀, 0.0], domain_symbol=:Γs)],
        resonator_domains = [G.ResonatorDomain(location=[resonator_x[k], 0.0], delta_symbol=Symbol("δ$k")) for k in 1:Nᵣ_ACTIVE],
    )

    potential = P.PotentialFlow(
      g=g,
      ρw=ρw,
      sea_state=inc_wave.sea_state,
      boundary_conditions=[
        P.RadiationBC(domain=:dΓin),
        P.RadiationBC(domain=:dΓout),
        P.PrescribedInletPotentialBC(domain=:dΓin, forcing=f_in, quantity=:traction),
      ],
      fe=PH.FESpaceConfig(
        order=2,
        vector_type=Vector{ComplexF64}
      ),
      space_domain_symbol=:Ω
    )
    free_surface = P.FreeSurface(
      ρw=ρw,
      g=g,
      βₕ=0.5,
      fe=PH.FESpaceConfig(order=2, vector_type=Vector{ComplexF64}),
      space_domain_symbol=:Γκ,
    )
    beam = P.EulerBernoulliBeam(
      L=L,
      mᵨ=mᵨ,
      EIᵨ=EIᵨ,
      g=g,
      symbol=:w,
      fe=PH.FESpaceConfig(
        order=2,
        vector_type=Vector{ComplexF64},
        γ=6.0
      ),
      space_domain_symbol=:Γs
    )

    # Density-graded population: Mᵣ and Kᵣ are BOTH constant across the
    # active set (unlike the frequency-graded script, which varies Kᵣ per
    # cell). If your installed resonator_array requires all three
    # positional arguments to be same-length vectors, this already is one.
    resonators = P.resonator_array(Nᵣ_ACTIVE, fill(Mᵣ, Nᵣ_ACTIVE), fill(Kᵣ, Nᵣ_ACTIVE), zeros(Nᵣ_ACTIVE); delta_domain_symbols=[Symbol("δ$k") for k in 1:Nᵣ_ACTIVE], ρw=ρw, fe=PH.FESpaceConfig(space_type=P.VectorValue{1,Float64}, vector_type=Vector{ComplexF64}))
    entities = P.PhysicsParameters[potential, free_surface, beam, resonators]

    config = PH.FreqDomainConfig(ω=ωinc)
    problem = S.build_problem(tank, entities, config)
    println(" Problem built, about to simulate case ωinc=$ωinc.")
    ϕₕ,κₕ,ηₕ,ri = S.simulate(problem).solution
    trians = S.get_triangulations(problem)

    # --- Sample the plate deflection along the plate (local coordinate
    # ξ = x - x₀), for the elevation-vs-frequency map built in the sweep
    # below. Identical to the frequency-graded script.
    ξs = range(1e-6, Ls - 1e-6, length=NX_FIELD_SAMPLE)
    pts = [Point(x₀ + ξ, 0.0) for ξ in ξs]
    println(" Evaluating solution for case ωinc=$ωinc.")
    ηvals = ηₕ(pts)  # use evaluate(ηₕ, pts) if ηₕ doesn't broadcast

    if vtk_output
      writevtk(trians[:Ω], "tmp_vol_density_graded_$ωinc.vtu", cellfields=["ϕ" => real(ϕₕ),"ϕin"=>x->real(inc_wave.ϕin(x))])
      writevtk(trians[:Γκ], "tmp_surf_density_graded_$ωinc.vtu", cellfields=["κ" => abs(κₕ),"ηin" => x->abs(inc_wave.ηin(x))],nsubcells=4)
      writevtk(trians[:Γη], "tmp_struct_density_graded_$ωinc.vtu", cellfields=["η" => abs(ηₕ)],nsubcells=4)
    end

    xs_in  = PP.probe_positions(x₀, :upwave,   λ, H_margin; margin=PROBE_MARGIN, n=N_PROBES)
    xs_out = PP.probe_positions(x₀+Ls, :downwave, λ, H_margin; margin=PROBE_MARGIN, n=N_PROBES)
    @assert minimum(xs_in) >= 0 && maximum(xs_out) <= L "Probe positions fell outside the tank (λ=$λ, L=$L) — sizing invariant broken, please report."
    rt = PP.reflection_transmission_coefficients(κₕ, xs_in, xs_out, inc_wave.k; y=0, η0=η0)
    println("Results for ωinc = $ωinc rad/s, λ = $λ m:")
    println(rt)

    rt_nt = NamedTuple{fieldnames(typeof(rt))}(getfield.(Ref(rt), fieldnames(typeof(rt))))
    return merge(rt_nt, (; ξs=collect(ξs), ηvals=ηvals))
  end
end # @everywhere

# --- everything below runs only on the master process ---
using Plots

# Warm-up on every worker (JIT-compile latency paid once, up front).
@everywhere main(1, 10.0, false)

# Sweep over the same range as the frequency-graded script, for direct
# comparison (a single fixed ωᵣ means the interesting structure is
# concentrated near ω ≈ ω_r rather than spread over a band — narrow this
# range around ω_r if you only care about the local resonance).
ωinc_values = 1.0:0.01:12.0
results = pmap(ωinc -> main(30, ωinc, false), ωinc_values)

Rs = [r.R for r in results]
Ts = [r.T for r in results]
field_matrix = permutedims(reduce(hcat, (abs.(r.ηvals) for r in results)))  # rows=ω, cols=ξ
ξs_ref = range(1e-6, Ls - 1e-6, length=NX_FIELD_SAMPLE)

# --- Plot 1: R and T (log scale) ---
plt = plot(ωinc_values, abs.(Rs), label="R", yscale=:log10)
plot!(plt, ωinc_values, abs.(Ts), label="T", yscale=:log10)
savefig(plt, "R_T_coefficients_density_graded.png")

# --- Plot 2: R^2+T^2 energy-conservation check ---
plt2 = plot(ωinc_values, [abs(Rs[i])^2 + abs(Ts[i])^2 for i in 1:length(Rs)], label="R^2+T^2")
savefig(plt2, "R2_T2_coefficients_density_graded.png")

# --- Plot 3: |R| on a linear scale. A single vertical line marks ω_r
# (unlike the frequency-graded case, there is no separated band — the
# resonance sits at one fixed frequency and grading only controls how
# strongly it is expressed / how it decays spatially). ---
plt3 = plot(ωinc_values, abs.(Rs), label="|R|", ylim=(0,1.05),
            xlabel="ω (rad/s)", ylabel="|R|",
            title="Reflection amplitude, density-graded array")
vline!(plt3, [ω_r], linestyle=:dash, color=:black, label="ωᵣ = $(ω_r) rad/s")
savefig(plt3, "R_amplitude_density_graded.png")

# --- Plot 4: plate elevation |η(ξ)| vs (ξ, ω). The reference line is now
# HORIZONTAL at ω = ω_r (the resonance frequency does not vary spatially
# under density grading; only the coupling strength does). ---
plt4 = heatmap(ξs_ref, ωinc_values, field_matrix,
               xlabel="ξ = x - x₀ (m), distance along plate",
               ylabel="ω (rad/s)",
               title="Plate elevation |η| vs frequency (density-graded array)",
               color=:viridis)
hline!(plt4, [ω_r], linestyle=:dash, linewidth=2, color=:white, label="ωᵣ (fixed)")
savefig(plt4, "eta_field_vs_freq_density_graded.png")

# --- Plot 5: diagnostic — realised resonator occupation vs the target
# continuum density n(ξ)/n₀ = e^{-ξ/L_grade}. Lets you visually confirm
# the thinning pattern (which candidate cells were kept vs skipped)
# tracks the intended exponential law. ---
ζ_active  = [Lᵣ/2 + (i-1)*Lᵣ for i in ACTIVE_IDXS]
ζ_removed = [Lᵣ/2 + (i-1)*Lᵣ for i in setdiff(1:Nᵣ_MAX, ACTIVE_IDXS)]
ξ_curve = range(0, Ls, length=400)
density_curve = exp.(-ξ_curve ./ L_grade)
plt5 = plot(ξ_curve, density_curve, label="n(ξ)/n₀ (target)",
            xlabel="ξ = x - x₀ (m)", ylabel="normalized density / occupancy",
            ylim=(-0.05, 1.05),
            title="Resonator population thinning (density grading)")
scatter!(plt5, ζ_active, zeros(length(ζ_active)), label="occupied cell", marker=:vline, markersize=8, color=:black)
scatter!(plt5, ζ_removed, zeros(length(ζ_removed)), label="thinned-out cell", marker=:x, markersize=4, color=:red)
savefig(plt5, "resonator_occupation_density_graded.png")

println("Density-graded run complete: $(Nᵣ_ACTIVE) of $(Nᵣ_MAX) candidate cells populated.")