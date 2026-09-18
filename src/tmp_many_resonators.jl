module TMP_MANY_RESONATORS

using HydroElasticFEM
using Gridap
import WaveSpec as WS
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.PostProcessing as PP
using HydroElasticFEM: map_vertical_GP_for_const_dep
using Plots

Ls = 30.0
H = 10.0
g = WS.PhysicalConstants.g#g = 9.8
mᵨ = 0.01
EIᵨ = 0.05*g
ωᵣ = 10.0
Mᵣ = 10.0
Kᵣ = Mᵣ * ωᵣ^2
ρw = 1000.0
nᵣ = 30
Lᵣ = 1.0 # Unit cell length for resonator array
# ωinc = 5.0
η0 = 0.01

# Probe layout for R/T post-processing — kept in sync between the
# domain-sizing block and the PP.probe_positions() calls below, so the
# probes are guaranteed (by construction, see main()) to land inside the
# tank at every frequency in the sweep, not just the ones that happened
# to be tried.
const N_PROBES     = 4
const PROBE_MARGIN = 1.5   # × H standoff from the structure edge (evanescent decay)

function incident_wave(ω::Real, η0::Real)
  k = WS.AiryWaves.solve_wavenumber(ω, H)
  T = 2π / ω
  spec = WS.ContinuousSpectrums.RegularWave(2*η0, T)
  ds = WS.SpectralSpreading.DiscreteSpectralSpreading(spec; mess=false)
  spread = WS.AngularSpreading.DiscreteAngularSpreading(0.0)
  sea_state = WS.AiryWaves.AiryState(ds, spread, 1, 1, [ω], [k], [0.0], H, 1)
  # z is SWL-referenced here (z=0 at the free surface), consistent with
  # domain_map()'s shifted mesh. The textbook form is
  #   ϕ = -i(g η0/ω) cosh(k(z+H))/cosh(kH).
  # Rewriting the prefactor as η0*ω/k (via ω² = g k tanh(kH)) turns
  # cosh(kH) into sinh(kH) in the denominator — using cosh(kH) with the
  # ω/k prefactor (as before) silently scales the whole incident forcing
  # by an extra factor of tanh(kH): ≈1 (harmless) for kH ≳ 3, but ~16%
  # low at kH≈1.2 (ω≈1 rad/s here) — exactly the low end of [1,12].
  # The vertical velocity also needs sinh(k(z+H)), not cosh(k(z+H)): as
  # written before it gave nonzero flow through the seabed at z=-H. It
  # happens to have no effect on this script (f_in's dot product with
  # the inlet normal (-1,0) only ever touches the horizontal component)
  # but is fixed here since vin is otherwise mislabelled for any future
  # use (diagnostics, a differently-oriented boundary, ...).
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

function _damping(k)
    μ0 = 2.5
    μ1_in(x) = μ0 * (1.0 - sin(pi / 2 * (x[1]) / Ld)) * (x[1] <= Ld)
    μ1_out(x) = μ0 * (1.0 - cos(pi / 2 * (x[1] - (L - Ld)) / Ld)) * (x[1] >= (L - Ld))
    μ2_in(x) = μ1_in(x) * k
    μ2_out(x) = μ1_out(x) * k
    return (; μ1_in, μ1_out, μ2_in, μ2_out)
end

function main(n::Int,ωinc::Real,vtk_output::Bool=false)

  inc_wave = incident_wave(ωinc, η0)
  damp = _damping(inc_wave.k)
  f_in(x) = (inc_wave.vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * inc_wave.k * inc_wave.ϕin(x)

  λ = 2π / inc_wave.k
  println("Running simulation for ωinc = $ωinc rad/s, k = $(inc_wave.k) rad/m, λ = $λ m")

  # Domain size and discretization depend on the wavelength:
  # - 20 elements per wavelength in x-direction
  # - n elements in y-direction
  # - Structure Ls m long (fixed, independent of λ) — the margin on each
  #   side must clear (a) evanescent-mode decay near the structure /
  #   resonators, set by the depth H (frequency-independent — this is
  #   what dominates at HIGH ω, where a fixed "N wavelengths" margin
  #   shrinks to nothing and used to leave the structure sticking out of
  #   the tank), and (b) the probe spread PP.probe_positions itself
  #   needs, which does scale with λ (dominates at LOW ω). Reusing
  #   PP.suggest_probe_offsets for the sizing — rather than an
  #   independent "3 wavelengths" rule — keeps this block and the probe
  #   calls below exactly in sync.
  standoff = PROBE_MARGIN * H + maximum(PP.suggest_probe_offsets(λ; n=N_PROBES))
  x₀ = 1.25 * standoff   # 25% buffer so probes aren't flush with dΓin/dΓout
  L  = 2 * x₀ + Ls

  h  = min(λ / 20, Lᵣ / 2)
  nx = ceil(Int, L / h)

  tank = G.TankDomain(
      L = L, H = H, nx = nx, ny = n,
      map=domain_map(n),
      is_periodic = (false, false),
      # damping_zones = [
      #     G.DampingZone(L = Ld, x₀ = [0.0, 0.0], domain_symbol = :Γ_d_in),
      #     G.DampingZone(L = Ld, x₀ = [L - Ld, 0.0], domain_symbol = :Γ_d_out),
      # ],
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
      # P.DampingZoneBC(
      #     domain = :dΓ_d_in,
      #     μ₁ = damp.μ1_in,
      #     μ₂ = damp.μ2_in,
      #     η_in = inc_wave.ηin,
      #     vz_in = x -> inc_wave.vin(x)⋅VectorValue(0.0, 1.0),
      # ),
      # P.DampingZoneBC(
      #     domain = :dΓ_d_out,
      #     μ₁ = damp.μ1_out,
      #     μ₂ = damp.μ2_out,
      #     η_in = (x -> 0.0 + 0.0im),
      #     vz_in = (x -> 0.0 + 0.0im),
      # ),
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

  resonators = P.resonator_array(nᵣ, Mᵣ, Kᵣ, 0.0; delta_domain_symbols=[Symbol("δ$i") for i in 1:nᵣ], ρw=ρw, fe=PH.FESpaceConfig(space_type=P.VectorValue{1,Float64}, vector_type=Vector{ComplexF64}))
  entities = P.PhysicsParameters[potential, free_surface, beam, resonators]
  
  config = PH.FreqDomainConfig(ω=ωinc)
  problem = S.build_problem(tank, entities, config)
  ϕₕ,κₕ,ηₕ,ri = S.simulate(problem).solution
  trians = S.get_triangulations(problem)
  if vtk_output
    writevtk(trians[:Ω], "tmp_vol.vtu", cellfields=["ϕ" => real(ϕₕ),"ϕin"=>x->real(inc_wave.ϕin(x))],nsubcells=4)
    writevtk(trians[:Γκ], "tmp_surf.vtu", cellfields=["κ" => abs(κₕ),"ηin" => x->abs(inc_wave.ηin(x)), "μ1in" => damp.μ1_in, "μ1out" => damp.μ1_out],nsubcells=4)
    writevtk(trians[:Γη], "tmp_struct.vtu", cellfields=["η" => abs(ηₕ)],nsubcells=4)
  end

  # Post-processing reflected/transmitted wave amplitudes.
  # Same N_PROBES/PROBE_MARGIN as the domain-sizing block above — by
  # construction (see there) these always land inside [0, L]; the
  # assertion is a cheap safety net in case that invariant is ever
  # broken by a future change to either block.
  xs_in  = PP.probe_positions(x₀, :upwave,   λ, H; margin=PROBE_MARGIN, n=N_PROBES)
  xs_out = PP.probe_positions(x₀+Ls, :downwave, λ, H; margin=PROBE_MARGIN, n=N_PROBES)
  @assert minimum(xs_in) >= 0 && maximum(xs_out) <= L "Probe positions fell outside the tank (λ=$λ, L=$L) — sizing invariant broken, please report."
  rt = PP.reflection_transmission_coefficients(κₕ, xs_in, xs_out, inc_wave.k; y=0, η0=η0)
  println("Results for ωinc = $ωinc rad/s, λ = $λ m:")
  println(rt)
  
  return rt
end

# main(4)
main(10, 10.0, false)  # Warm-up example call with ωinc = 5.0 rad/s

# Sweep across frequencies range 1 to 12 rad/s
Rs = ComplexF64[]
Ts = ComplexF64[]
ωinc_values = 1.0:0.1:12.0
for ωinc in ωinc_values
  rt = main(10, ωinc, false)
  push!(Rs, rt.R)
  push!(Ts, rt.T)
end

# Plot Transmission and Reflection coefficients (log scale) in single plot
plt = plot(ωinc_values, abs.(Rs), label="R", yscale=:log10)
plot!(plt, ωinc_values, abs.(Ts), label="T", yscale=:log10)
savefig(plt, "R_T_coefficients.png")

# Plot R^2+T^2 to check energy conservation
plt2 = plot(ωinc_values, [abs(Rs[i])^2 + abs(Ts[i])^2 for i in 1:length(Rs)], label="R^2+T^2")
savefig(plt2, "R2_T2_coefficients.png")

end