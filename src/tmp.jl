module TMP

using HydroElasticFEM
using Gridap
import WaveSpec as WS
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S
import HydroElasticFEM.ParameterHandler as PH
using HydroElasticFEM: map_vertical_GP_for_const_dep

L = 130.0
Ld = 30.0
Ls = 30.0
x₀ = 50.0
H = 10.0
g = WS.PhysicalConstants.g#9.80665
mᵨ = 0.01
EIᵨ = 0.05*g
ωᵣ = 5.0
Mᵣ = 10.0
Kᵣ = Mᵣ * ωᵣ^2
ρw = 1025.0

function incident_wave(; H0::Real, ω::Real, η0::Real, α::Real)
  T = 2π / ω
  spec = WS.ContinuousSpectrums.RegularWave(2*η0, T)
  ds = WS.SpectralSpreading.DiscreteSpectralSpreading(spec; mess=false)
  spread = WS.AngularSpreading.DiscreteAngularSpreading(α)
  k = WS.AiryWaves.solve_wavenumber(ω, H)
  println(k*g*tanh(k*H) - ω^2)
  sea_state = WS.AiryWaves.AiryState(ds, spread, 1, 1, [ω], [k], [α], H, 1)
  wave(x) = WS.AiryWaves.generate_sea(sea_state, [x[1]], [0.0], [x[2]], [0.0], vars=[:η, :ϕ, :u, :w])
  ηin(x) = η0*exp(im*k*x[1])#wave(x)[:η][1]
  # ϕin(x) = wave(x)[:ϕ][1]
  ϕin(x) = -im*(η0*ω/k)*(cosh(k*(x[2]+H)) / cosh(k*H))*exp(im*k*x[1])
  # vin(x) = VectorValue(wave(x)[:u][1], wave(x)[:w][1])
  vin(x) = VectorValue((η0*ω)*(cosh(k*(x[2]+H)) / cosh(k*H))*exp(im*k*x[1]),
                       -im*(η0*ω)*(cosh(k*(x[2]+H)) / cosh(k*H))*exp(im*k*x[1]))
  return (; sea_state, ηin, ϕin, vin)
end

function domain_map(ny)
  x -> VectorValue( x[1],
    map_vertical_GP_for_const_dep(x[2] - H, 1.08, ny, H; dbgmsg=false),
  )
end

function _damping(k)
    μ0 = 0.1
    μ1_in(x) = μ0 * (1.0 - sin(pi / 2 * (x[1]) / Ld)) * (x[1] <= Ld)
    μ1_out(x) = μ0 * (1.0 - cos(pi / 2 * (x[1] - (L - Ld)) / Ld)) * (x[1] >= (L - Ld))
    μ2_in(x) = μ1_in(x) * k
    μ2_out(x) = μ1_out(x) * k
    return (; μ1_in, μ1_out, μ2_in, μ2_out)
end

function main(n)

  inc_wave = incident_wave(H0=H, ω=ωᵣ, η0=0.01, α=0.0)
  damp = _damping(inc_wave.sea_state.k[1])
  f_in(x) = (inc_wave.vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * inc_wave.sea_state.k[1] * inc_wave.ϕin(x)

  tank = G.TankDomain(
      L = L, H = H, nx = 13*n, ny = n,
      map=domain_map(n),
      is_periodic = (false, false),
      # damping_zones = [
      #     G.DampingZone(L = Ld, x₀ = [0.0, 0.0], domain_symbol = :Γ_d_in),
      #     # G.DampingZone(L = Ld, x₀ = [L - Ld, 0.0], domain_symbol = :Γ_d_out),
      # ],
      structure_domains = [G.StructureDomain(L=Ls, x₀=[x₀, 0.0], domain_symbol=:Γs)],
      resonator_domains = [G.ResonatorDomain(location=[x₀+10.0,0.0],delta_symbol=:δ1), G.ResonatorDomain(location=[x₀+20.0,0.0],delta_symbol=:δ2)],
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

  resonators = P.resonator_array(2, Mᵣ, Kᵣ, 0.0; delta_domain_symbols=[:δ1, :δ2], ρw=ρw, fe=PH.FESpaceConfig(space_type=P.VectorValue{1,Float64}, vector_type=Vector{ComplexF64}))
  entities = P.PhysicsParameters[potential, free_surface, beam, resonators]
  # entities = P.PhysicsParameters[potential, free_surface]

  config = PH.FreqDomainConfig(ω=ωᵣ)
  problem = S.build_problem(tank, entities, config)
  ϕₕ,κₕ,ηₕ,r1,r2 = S.simulate(problem).solution
  # ϕₕ,κₕ = S.simulate(problem).solution
  trians = S.get_triangulations(problem)
  writevtk(trians[:Ω], "tmp_vol.vtu", cellfields=["ϕ" => real(ϕₕ),"ϕin"=>x->real(inc_wave.ϕin(x))],nsubcells=4)
  writevtk(trians[:Γκ], "tmp_surf.vtu", cellfields=["κ" => abs(κₕ),"ηin" => x->abs(inc_wave.ηin(x)), "μ1in" => damp.μ1_in, "μ1out" => damp.μ1_out],nsubcells=4)
  writevtk(trians[:Γη], "tmp_struct.vtu", cellfields=["η" => abs(ηₕ), "r1" => real(r1), "r2" => real(r2)],nsubcells=4)
end

main(4)
main(40)

end

