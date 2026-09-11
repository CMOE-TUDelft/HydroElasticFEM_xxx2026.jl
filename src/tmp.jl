using HydroElasticFEM
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S
import HydroElasticFEM.ParameterHandler as PH

L = 4.0
H = 1.0
n = 4

stop

tank = G.TankDomain(
    L = L, H = H, nx = 2*n, ny = n,
    is_periodic = (true, false),
    structure_domains = [G.StructureDomain(L=L, x₀=[0.0, H], domain_symbol=:Γs)],
)

model = G.build_model(tank)
trians = G.build_triangulations(tank, model)
println("trian keys: ", keys(trians.data))

potential = P.PotentialFlow(g=9.81, fe=PH.FESpaceConfig(order=2, vector_type=Vector{Float64}), space_domain_symbol=:Ω)
beam = P.EulerBernoulliBeam(L=L, mᵨ=0.01, EIᵨ=0.49, g=9.81, symbol=:w,
    fe=PH.FESpaceConfig(order=2, vector_type=Vector{Float64}, γ=6.0), space_domain_symbol=:Γs)

positions = [P.VectorValue(0.5, H), P.VectorValue(2.5, H)]
M_r = 10.0 * 1000.0
omega_r = 10.0
K_r = M_r * omega_r^2
resn0 = P.resonator_array(2, M_r, K_r, 0.0, positions; ρw=1000.0)
resn = [P.ResonatorSingle(M=r.M, K=r.K, C=r.C, ρw=r.ρw, XZ=r.XZ,
        fe=PH.FESpaceConfig(space_type=P.VectorValue{1,Float64}, vector_type=Vector{Float64})) for r in resn0]

entities = Any[potential, beam, resn]
X, Y, fmap = S.FESpaceAssembly.build_fe_spaces(entities, trians, PH.TimeDomainConfig())
println("fmap: ", fmap)

deg = S.get_integration_degrees(trians, P.PhysicsParameters[potential, beam])
dom = G.get_integration_domains(trians; degree=deg)
Γtop = trians[:Γη]
dom[:δ_p] = [DiracDelta(Γtop, [Gridap.Point(x[1], x[2])]) for x in positions]
println("dom ready, keys: ", keys(dom))
