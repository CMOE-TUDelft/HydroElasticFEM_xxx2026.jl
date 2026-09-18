# Frequency-domain response of a floating metaplate with GRADED resonator
# frequencies, following Liu, Farhat, Bagci, Guenneau & Wu (2025), J. Fluid
# Mech. 1020, A50, §5.1-5.2 ("Graded floating metaplate: broadband wave
# reflector with customisable working frequency range").
#
# This is a direct adaptation of the attached uniform-array script
# (TMP_MANY_RESONATORS): same geometry, material and mesh conventions, only
# the resonator array is changed from a single (Mᵣ, Kᵣ) pair to a per-
# resonator (Mᵣ, Kᵣ(i)) array so that each resonator has its own natural
# frequency ωᵣ(i).
#
# PARALLELISED WITH Distributed.jl. Previously this file was wrapped in
# `module GRADED_RESONATORS_FREQ ... end`. That wrapper is dropped here:
# Distributed workers are separate OS processes with no shared memory, so
# `main` and everything it needs (package imports, constants,
# incident_wave, domain_map) have to be pushed to them explicitly via
# `@everywhere`, and `pmap` then looks them up as plain top-level (Main)
# names on each worker. Keeping a custom module would mean replicating
# that same module on every worker just to match the master's
# `GRADED_RESONATORS_FREQ.main` — unnecessary complexity for a driver
# script (as opposed to a package).

using Distributed

# Number of worker processes. Defaults to "all cores but one" (leaving one
# free for the OS / master process); tune to your machine, or to fewer if
# your meshes are large enough that a single solve already saturates
# several cores internally. Guarded so re-`include`-ing this file in the
# same session doesn't keep adding more workers on top of existing ones.
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

  # Each worker is its own OS process; without this, BLAS inside each
  # worker's FEM solve would itself try to spawn multiple threads,
  # oversubscribing cores across workers.
  BLAS.set_num_threads(1)

  Ls = 30.0
  H = 10.0
  g = WS.PhysicalConstants.g
  mᵨ = 0.01
  EIᵨ = 0.05*g
  Mᵣ = 10.0            # resonator mass, kept CONSTANT across the array (§5.1)
  ρw = 1000.0
  nᵣ = 30
  Lᵣ = 1.0
  η0 = 0.01

  # ---------------------------------------------------------------------
  # Graded resonant frequency, §5.1 / figure 8(a) / Table 1.
  #
  # The paper grades ωᵣ linearly from 8 to 11 rad/s across the N=30 unit
  # cells, with the mass held fixed and the spring stiffness varied instead:
  #     sᵣ(i) = Mᵣ · ωᵣ(i)²                                          (5.1)
  #
  # ORIENTATION: in the paper's own figure the incident wave arrives from
  # the +x (downwave) side and ωᵣ DEcreases with x (11 rad/s at the left,
  # 8 rad/s at the right; fig. 8a). The tank re-used from the attached
  # script instead launches the incident wave from the x=0 (dΓin/upwave)
  # side. To reproduce the SAME rainbow-trapping physics (the wave starts
  # in a region where it is in the passband, ωᵣ(x) < ω_inc, and is
  # progressively trapped as ωᵣ(x) rises to meet ω_inc) the grading below
  # is mirrored: ωᵣ increases with x, i.e. LOW near the inlet and HIGH
  # near the outlet. Swap ω_res_lo/ω_res_hi if your incident wave enters
  # from the opposite side.
  # ---------------------------------------------------------------------
  const ω_res_lo = 8.0    # rad/s, resonator nearest dΓin  (i = 1)
  const ω_res_hi = 11.0   # rad/s, resonator nearest dΓout (i = nᵣ)
  ωᵣ_vec = collect(range(ω_res_lo, ω_res_hi, length=nᵣ))
  Kᵣ_vec = Mᵣ .* ωᵣ_vec .^ 2

  # Probe layout for R/T post-processing (unchanged from the attached script)
  const N_PROBES     = 4
  const PROBE_MARGIN = 1.5

  # Number of points sampled along the plate for the fig.10-style
  # elevation-vs-frequency map (see the field-sampling block in main()
  # and plot 4 at the bottom of the sweep).
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
    x -> VectorValue( x[1],
      map_vertical_GP_for_const_dep(x[2] - H, 1.08, ny, H; dbgmsg=false),
    )
  end

  function main(n::Int, ωinc::Real, vtk_output::Bool=false)

    inc_wave = incident_wave(ωinc, η0)
    f_in(x) = (inc_wave.vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * inc_wave.k * inc_wave.ϕin(x)

    λ = 2π / inc_wave.k
    println("Graded array — ωinc = $ωinc rad/s, k = $(inc_wave.k) rad/m, λ = $λ m")

    # Depth used for the evanescent-decay margin. PROBE_MARGIN*H alone is
    # frequency-INDEPENDENT (sized for the worst case: low ω, slow decay),
    # so at high ω the domain length L stayed ~constant while h = λ/20
    # kept shrinking, driving nx = L/h up. Physically the dominant
    # evanescent wavenumber grows with ω, so the true decay length
    # shrinks at high frequency — capping the depth at min(H, λ) lets the
    # margin (and hence L) shrink accordingly at high ω, while leaving the
    # low-frequency (λ ≳ H) sizing exactly as before. H_margin is reused
    # below for the actual probe placement so the two stay in sync (same
    # reasoning as the original "keeps this block and the probe calls
    # below exactly in sync" comment).
    H_margin = min(H, λ)
    standoff = ceil(PROBE_MARGIN * H_margin + maximum(PP.suggest_probe_offsets(λ; n=N_PROBES)))
    h  = min(λ / 20, Lᵣ / 2)
    nxs = ceil(Int, (Lᵣ/2) / h)*ceil(Int, Ls/(Lᵣ/2))
    hs = Ls / nxs
    standoff_ = ceil(Int, standoff/hs)*hs
    L = 2*standoff_ + Ls
    x₀ = standoff_
    nx = Int(L / hs)
    println("Tank length L = $L m, nx = $nx, ny = $n, h = $h m")

    tank = G.TankDomain(
        L = L, H = H, nx = nx, ny = n,
        map=domain_map(n),
        is_periodic = (false, false),
        structure_domains = [G.StructureDomain(L=Ls, x₀=[x₀, 0.0], domain_symbol=:Γs)],
        resonator_domains = [G.ResonatorDomain(location=[x₀+Lᵣ/2+(i-1)*Lᵣ,0.0],delta_symbol=Symbol("δ$i")) for i in 1:nᵣ],
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

    # Graded array: pass the per-resonator stiffness vector Kᵣ_vec instead of
    # a scalar Kᵣ. Mass Mᵣ stays a scalar since it is uniform across the
    # array (§5.1). If your installed resonator_array requires all three
    # positional arguments to be same-length vectors, use instead:
    #   P.resonator_array(nᵣ, fill(Mᵣ, nᵣ), Kᵣ_vec, zeros(nᵣ); ...)
    resonators = P.resonator_array(nᵣ, fill(Mᵣ, nᵣ), Kᵣ_vec, zeros(nᵣ); delta_domain_symbols=[Symbol("δ$i") for i in 1:nᵣ], ρw=ρw, fe=PH.FESpaceConfig(space_type=P.VectorValue{1,Float64}, vector_type=Vector{ComplexF64}))
    entities = P.PhysicsParameters[potential, free_surface, beam, resonators]

    config = PH.FreqDomainConfig(ω=ωinc)
    problem = S.build_problem(tank, entities, config)
    println(" Problem built, about to simulate case ωinc=$ωinc.")
    ϕₕ,κₕ,ηₕ,ri = S.simulate(problem).solution
    trians = S.get_triangulations(problem)

    # --- Sample the plate deflection along the plate, for the fig.10-style
    # elevation-vs-frequency map built in the sweep below. We sample in the
    # LOCAL plate coordinate ξ = x - x₀ ∈ [0, Ls] rather than absolute x,
    # because x₀ itself shifts with frequency (the tank is resized per
    # ωinc above) — ξ is what stays comparable across the whole sweep.
    # Points are nudged 1e-6 m in from each edge to avoid landing exactly
    # on the Γs/open-water boundary. If ηₕ.(pts) doesn't broadcast on your
    # Gridap version, use `evaluate(ηₕ, pts)` instead.
    ξs = range(1e-6, Ls - 1e-6, length=NX_FIELD_SAMPLE)
    pts = [Point(x₀ + ξ, 0.0) for ξ in ξs]
    # With h=0.5 m and NX_FIELD_SAMPLE points spread over Ls=30 m (~0.15 m
    # spacing), many sample points land close to *some* internal element
    # boundary, not just the two plate edges — the default 1-nearest-
    # vertex KD-tree search can then fail to locate the containing cell.
    # Wrapping ηₕ as an Interpolable with a larger search neighbourhood
    # fixes this (this is exactly what Gridap's own error message
    # suggests). Bump num_nearest_vertices further if it still happens.
    # ηₕ_interp = Gridap.CellData.Interpolable(ηₕ; searchmethod=Gridap.CellData.KDTreeSearch(num_nearest_vertices=3))
    # ηvals = ηₕ_interp.(pts)
    println(" Evaluating solution for case ωinc=$ωinc.")
    ηvals = ηₕ(pts)  # use evaluate(ηₕ, pts) if ηₕ doesn't broadcast

    if vtk_output
      writevtk(trians[:Ω], "tmp_vol_graded_$ωinc.vtu", cellfields=["ϕ" => real(ϕₕ),"ϕin"=>x->real(inc_wave.ϕin(x))])
      writevtk(trians[:Γκ], "tmp_surf_graded_$ωinc.vtu", cellfields=["κ" => abs(κₕ),"ηin" => x->abs(inc_wave.ηin(x))],nsubcells=4)
      writevtk(trians[:Γη], "tmp_struct_graded_$ωinc.vtu", cellfields=["η" => abs(ηₕ)],nsubcells=4)
    end

    xs_in  = PP.probe_positions(x₀, :upwave,   λ, H_margin; margin=PROBE_MARGIN, n=N_PROBES)
    xs_out = PP.probe_positions(x₀+Ls, :downwave, λ, H_margin; margin=PROBE_MARGIN, n=N_PROBES)
    @assert minimum(xs_in) >= 0 && maximum(xs_out) <= L "Probe positions fell outside the tank (λ=$λ, L=$L) — sizing invariant broken, please report."
    rt = PP.reflection_transmission_coefficients(κₕ, xs_in, xs_out, inc_wave.k; y=0, η0=η0)
    println("Results for ωinc = $ωinc rad/s, λ = $λ m:")
    println(rt)

    # merge keeps rt.R / rt.T working exactly as before at the call site,
    # while also exposing the along-plate field sampled above. rt is a
    # struct (ReflectionTransmissionResult), not a NamedTuple, so it has
    # to be converted first — done generically via fieldnames/getfield so
    # this doesn't hard-code which fields the struct carries.
    rt_nt = NamedTuple{fieldnames(typeof(rt))}(getfield.(Ref(rt), fieldnames(typeof(rt))))
    return merge(rt_nt, (; ξs=collect(ξs), ηvals=ηvals))
  end
end # @everywhere

# --- everything below runs only on the master process ---
using Plots   # only needed for the final plots, not on the workers

# Warm-up: compiled on EVERY worker (not just the master), so the first
# real job pmap sends to each one isn't also paying JIT-compile latency.
@everywhere main(1, 10.0, false)

# Sweep across frequencies range 1 to 12 rad/s — same range as figure 9.
#
# pmap distributes the calls across the worker pool and hands results
# back in the SAME ORDER as ωinc_values (unlike @distributed's default
# reduction), and load-balances dynamically — useful here since solve
# cost varies a lot across the sweep (nx grows at both the low- and
# high-ω ends, see the domain-sizing block in main()).
ωinc_values = 1.0:0.01:12.0
results = pmap(ωinc -> main(30, ωinc, false), ωinc_values)

Rs = [r.R for r in results]
Ts = [r.T for r in results]
field_matrix = permutedims(reduce(hcat, (abs.(r.ηvals) for r in results)))  # rows=ω, cols=ξ
# ξ grid for the along-plate sampling — identical every call inside main()
# since it only depends on Ls and NX_FIELD_SAMPLE (both fixed constants).
ξs_ref = range(1e-6, Ls - 1e-6, length=NX_FIELD_SAMPLE)

# --- Plot 1: R and T (log scale), same style as the attached script ---
plt = plot(ωinc_values, abs.(Rs), label="R", yscale=:log10)
plot!(plt, ωinc_values, abs.(Ts), label="T", yscale=:log10)
savefig(plt, "R_T_coefficients_graded.png")

# --- Plot 2: R^2+T^2 energy-conservation check ---
plt2 = plot(ωinc_values, [abs(Rs[i])^2 + abs(Ts[i])^2 for i in 1:length(Rs)], label="R^2+T^2")
savefig(plt2, "R2_T2_coefficients_graded.png")

# --- Plot 3: |R| on a linear scale, matching figure 9 of Liu et al. (2025) ---
# The shaded band marks the graded resonant-frequency range [8,11] rad/s,
# within which the paper reports near-total reflection.
plt3 = plot(ωinc_values, abs.(Rs), label="|R|", ylim=(0,1.05),
            xlabel="ω (rad/s)", ylabel="|R|",
            title="Reflection amplitude, graded array (cf. fig. 9)")
vspan!(plt3, [ω_res_lo, ω_res_hi], alpha=0.15, label="ωᵣ range (8–11)")
savefig(plt3, "R_amplitude_graded_fig9.png")

# --- Plot 4: plate elevation |η(ξ)| vs (ξ, ω), matching figure 10 ---
# Rows = incident frequency, columns = distance along the plate (from the
# upwave edge, ξ=0, to the downwave edge, ξ=Ls). The dashed line traces
# the local grading ωᵣ(ξ) (§5.1, linear 8→11 rad/s); the paper's
# "rainbow" signature is the diagonal ridge of high |η| that tracks this
# line — each frequency is amplified/trapped right where ω_inc ≈ ωᵣ(ξ),
# before being reflected.
plt4 = heatmap(ξs_ref, ωinc_values, field_matrix,
               xlabel="ξ = x - x₀ (m), distance along plate",
               ylabel="ω (rad/s)",
               title="Plate elevation |η| vs frequency (cf. fig. 10)",
               color=:viridis)
ωᵣ_of_ξ(ξ) = ω_res_lo + (ω_res_hi - ω_res_lo) * (ξ / Ls)   # linear grading guide line
plot!(plt4, collect(ξs_ref), ωᵣ_of_ξ.(collect(ξs_ref)),
      linestyle=:dash, linewidth=2, color=:white, label="local ωᵣ(ξ)")
savefig(plt4, "eta_field_vs_freq_graded_fig10.png")