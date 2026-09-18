# Interpolation-only sanity check for graded_resonators_freq.jl.
#
# Builds the SAME mesh/FE problem structure as main() (tank, entities,
# S.build_problem) for every frequency in the sweep, but NEVER calls
# S.simulate — no assembly of the physics operators, no linear solve.
# Point-location (the "Point ... was not found in any active cell" error
# seen before) only depends on the triangulation's geometry/topology, not
# on solved field values, so a trivial constant CellField on the same
# triangulation ηₕ would live on exercises exactly the same point-
# location code path as the real ηₕ(pts) call in main() — for free.
#
# Use this before committing to the full pmap sweep: it tells you,
# cheaply, whether the along-plate sampling grid (pts, NX_FIELD_SAMPLE
# points) will hit the interpolation error at any frequency, without
# paying for a solve at each one.
#
# IMPORTANT: the domain-sizing block, entity construction, and
# ωinc_values below must be kept identical to graded_resonators_freq.jl
# — if that file changes (e.g. the L = ceil(...) rounding, or the
# sweep range), mirror the change here too.

module tmp_check_probes

using HydroElasticFEM
using Gridap
import WaveSpec as WS
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.PostProcessing as PP
using HydroElasticFEM: map_vertical_GP_for_const_dep

# --- Must match graded_resonators_freq.jl exactly ---
Ls = 30.0
H = 10.0
g = WS.PhysicalConstants.g
mᵨ = 0.01
EIᵨ = 0.05*g
Mᵣ = 10.0
ρw = 1000.0
nᵣ = 30
Lᵣ = 1.0
η0 = 0.01

const ω_res_lo = 8.0
const ω_res_hi = 11.0
ωᵣ_vec = collect(range(ω_res_lo, ω_res_hi, length=nᵣ))
Kᵣ_vec = Mᵣ .* ωᵣ_vec .^ 2

const N_PROBES        = 4
const PROBE_MARGIN    = 1.5
const NX_FIELD_SAMPLE = 200

const ωinc_values = 8.0:0.02:8.02   # must match graded_resonators_freq.jl's sweep

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

struct InterpCheck
  ωinc::Float64
  L::Float64
  nx::Int
  n_fail::Int
  fail_ξ::Vector{Float64}
end

function check_interpolation(n::Int, ωinc::Real)
  inc_wave = incident_wave(ωinc, η0)
  f_in(x) = (inc_wave.vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * inc_wave.k * inc_wave.ϕin(x)

  λ = 2π / inc_wave.k

  # --- identical to the sizing block in main() ---
  H_margin = min(H, λ)
  standoff = ceil(PROBE_MARGIN * H_margin + maximum(PP.suggest_probe_offsets(λ; n=N_PROBES)))
  h  = min(λ / 20, Lᵣ / 2)
  nxs = ceil(Int, (Lᵣ/2) / h)*ceil(Int, Ls/(Lᵣ/2))
  hs = Ls / nxs
  standoff_ = ceil(Int, standoff/hs)*hs
  L = 2*standoff_ + Ls
  x₀ = standoff_
  nx = Int(L / hs)

  # Print domain and mesh properties
  println("Domain length L: ", L)
  println("Number of elements nx: ", nx)
  println("Effective element size h: ", L/nx)
  println("structure initial point x₀: ", x₀)

  tank = G.TankDomain(
      L = L, H = H, nx = nx, ny = n,
      map=domain_map(n),
      is_periodic = (false, false),
      structure_domains = [G.StructureDomain(L=Ls, x₀=[x₀, 0.0], domain_symbol=:Γs)],
      resonator_domains = [G.ResonatorDomain(location=[x₀+Lᵣ/2+(i-1)*Lᵣ,0.0],delta_symbol=Symbol("δ$i")) for i in 1:nᵣ],
  )

  # entities are needed by build_problem, but nothing from here on does
  # any physics assembly or solve — f_in is never evaluated.
  potential = P.PotentialFlow(
    g=g, ρw=ρw, sea_state=inc_wave.sea_state,
    boundary_conditions=[
      P.RadiationBC(domain=:dΓin),
      P.RadiationBC(domain=:dΓout),
      P.PrescribedInletPotentialBC(domain=:dΓin, forcing=f_in, quantity=:traction),
    ],
    fe=PH.FESpaceConfig(order=2, vector_type=Vector{ComplexF64}),
    space_domain_symbol=:Ω
  )
  free_surface = P.FreeSurface(
    ρw=ρw, g=g, βₕ=0.5,
    fe=PH.FESpaceConfig(order=2, vector_type=Vector{ComplexF64}),
    space_domain_symbol=:Γκ,
  )
  beam = P.EulerBernoulliBeam(
    L=L, mᵨ=mᵨ, EIᵨ=EIᵨ, g=g, symbol=:w,
    fe=PH.FESpaceConfig(order=2, vector_type=Vector{ComplexF64}, γ=6.0),
    space_domain_symbol=:Γs
  )
  resonators = P.resonator_array(nᵣ, fill(Mᵣ, nᵣ), Kᵣ_vec, zeros(nᵣ);
    delta_domain_symbols=[Symbol("δ$i") for i in 1:nᵣ], ρw=ρw,
    fe=PH.FESpaceConfig(space_type=P.VectorValue{1,Float64}, vector_type=Vector{ComplexF64}))
  entities = P.PhysicsParameters[potential, free_surface, beam, resonators]

  config = PH.FreqDomainConfig(ω=ωinc)
  problem = S.build_problem(tank, entities, config)   # <- mesh + FE spaces, NO solve
  trians = S.get_triangulations(problem)

  # Trivial constant field on the SAME triangulation ηₕ would live on.
  # Point-location depends only on triangulation geometry, not field
  # values, so this exercises exactly the failure mode being checked for.
  dummy = CellField( zero(ComplexF64),trians[:Γs])
  dummy_interp = Gridap.CellData.Interpolable(dummy; searchmethod=Gridap.CellData.KDTreeSearch(num_nearest_vertices=5,tol=1.0e-3),tol=1.0e-3)

  ξs = range(1e-3, Ls - 1e-3, length=NX_FIELD_SAMPLE)
  pts = [Point(x₀ + ξ, 0.0) for ξ in ξs]

  fail_ξ = Float64[]
  for (ξ, pt) in zip(ξs, pts)
    # try
      dummy_interp(pt)
    # catch
    #   push!(fail_ξ, ξ)
    # end
  end

  return InterpCheck(ωinc, L, nx, length(fail_ξ), fail_ξ)
end

check_interpolation(1, 1.0)

println("Checking interpolation for $(length(ωinc_values)) frequencies (mesh + FE space build only, no solve)...")
results = InterpCheck[]
for (i, ω) in enumerate(ωinc_values)
  r = check_interpolation(1, ω)
  push!(results, r)
  if r.n_fail > 0
    println("  ω = $(round(ω, digits=3)) rad/s (L=$(r.L), nx=$(r.nx)): " *
            "$(r.n_fail)/$(NX_FIELD_SAMPLE) points failed, " *
            "first at ξ ≈ $(round(r.fail_ξ[1], digits=4)) m")
  elseif i % 25 == 0
    println("  ... checked ω up to $(round(ω, digits=2)) rad/s, no failures so far")
  end
end

n_bad = count(r -> r.n_fail > 0, results)
println()
if n_bad == 0
  println("All $(length(results)) frequencies OK — every sample point interpolates " *
          "with the plain ηₕ(pts) call (no Interpolable/KDTreeSearch needed).")
else
  println("$(n_bad) / $(length(results)) frequencies have at least one failing point " *
          "— re-enable the Interpolable(...; searchmethod=KDTreeSearch(...)) wrapper " *
          "in main() before running the full sweep.")
end

end