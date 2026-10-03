module TMP3D_RESONATOR_HL

# Same problem as tmp_3D_resonator.jl (floating Kirchhoff-Love plate carrying a
# graded array of point mass-spring-damper resonators in a short-crested
# JONSWAP sea), assembled with the high-level API of HydroElasticFEM.jl:
# `ResonatorArray` + `G.ResonatorDomain` + build_problem/simulate.
# Tank, damping zones and incident sea are as in tmp_3D_hl.jl.

import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as S
import WaveSpec as WS
using Gridap
using TimerOutputs

to = TimerOutput()

# ─────────────────────────────────────────────────────────────────────────────
# Physical parameters
# ─────────────────────────────────────────────────────────────────────────────
H  = 30.0
ρ  = 1000.0
g  = 9.81

# Plate (Yago-type VLFS material, 60 × 50 m)
Lp = 60.0
Bp = 50.0
E  = 11.9e9
ρb = 256.25
ν  = 0.13
hb = 2.0

# Resonator array: nrx × nry resonators at the centres of equal plate patches
#   mass_ratio  total resonator mass / plate mass (M equal for all resonators)
#   ζr          damping ratio (Cᵢ = 2 ζr M ωrᵢ)
nrx, nry   = 6, 5
mass_ratio = 0.1
ζr         = 0.02

# Sea state
Hs = 3.0
Tp = 10.0
θ̄  = 0.0
σθ = π / 6
θmin, θmax = -π / 4, π / 4
fmin_fac = 0.7
μ₀ = 2.5

# ωr rises linearly in x across the plate, centred on the spectral peak ωp.
ωp = 2π / Tp
ωr_lo = 0.85 * ωp
ωr_hi = 1.15 * ωp

βₕ = 0.5
γₜ = 0.5
βₜ = 0.25
order = 2

# ─────────────────────────────────────────────────────────────────────────────
# Discretisation presets (see tmp_3D.jl for the rationale)
# ─────────────────────────────────────────────────────────────────────────────
const SMOKE = (; Δ = 10.0, nz = 1, zbase = 2.0, fmax_fac = 2.0, nf = 5, nθ = 3,
                 Δt = 0.2, T = 0.4, Δt_out = 0.2, seed = 1234, vtk_output = false)
const DEV = (; Δ = 10.0, nz = 3, zbase = 2.0, fmax_fac = 2.0, nf = 40, nθ = 11,
               Δt = 0.2, T = 400.0, Δt_out = 0.4, seed = 1234, vtk_output = true)
const PROD = (; Δ = 5.0, nz = 5, zbase = 2.0, fmax_fac = 2.5, nf = 50, nθ = 11,
                Δt = 0.1, T = 1400.0, Δt_out = 0.5, seed = 1234, vtk_output = true)

# ─────────────────────────────────────────────────────────────────────────────
# Wave scales, tank layout and damping (as in tmp_3D.jl)
# ─────────────────────────────────────────────────────────────────────────────

function wave_scales(fmax_fac)
  fp = 1 / Tp
  fmin, fmax = fmin_fac * fp, fmax_fac * fp
  k(f) = WS.AiryWaves.solve_wavenumber(2π * f, H)
  kp = k(fp)
  (; fp, fmin, fmax, kp, λp = 2π / kp, λmax = 2π / k(fmin), λmin = 2π / k(fmax),
     Tmin = 1 / fmax)
end

_snap(x, Δ) = Δ * ceil(x / Δ - 1e-9)

function tank_layout(ws, Δ)
  isinteger(Lp / Δ) && isinteger(Bp / Δ) ||
    error("Plate dimensions ($Lp × $Bp m) must be multiples of Δ = $Δ m.")
  Ld_in  = _snap(ws.λmax, Δ)
  Ld_out = _snap(2 * ws.λp, Δ)
  Ld_lat = _snap(ws.λp, Δ)
  G_up   = _snap(ws.λp, Δ)
  G_down = _snap(ws.λp, Δ)
  G_lat  = _snap(ws.λp / 2, Δ)
  L  = Ld_in + G_up + Lp + G_down + Ld_out
  W  = Bp + 2 * (G_lat + Ld_lat)
  xp₀ = Ld_in + G_up
  xd_out = L - Ld_out
  (; L, W, Ld_in, Ld_out, Ld_lat, G_up, G_down, G_lat, xp₀, xd_out,
     nx = round(Int, L / Δ), ny = round(Int, W / Δ))
end

# Locations and per-resonator (M, K, C, ωr): one at the centre of each plate patch.
function resonator_array(lay)
  N  = nrx * nry
  ωcol = nrx == 1 ? [ωp] : collect(range(ωr_lo, ωr_hi, length = nrx))
  idx  = [(ix, jy) for jy in 1:nry for ix in 1:nrx]
  locs = [[lay.xp₀ + (ix - 0.5) * Lp / nrx, -Bp / 2 + (jy - 0.5) * Bp / nry, 0.0]
          for (ix, jy) in idx]
  ωr = [ωcol[ix] for (ix, _) in idx]
  M  = fill(mass_ratio * ρb * hb * Lp * Bp / N, N)
  K  = M .* ωr .^ 2
  Cr = 2 .* ζr .* M .* ωr
  (; N, M, K, C = Cr, ωr, locs, δsyms = [Symbol("δ_r$i") for i in 1:N])
end

# Damping zones are numbered :dΓd_1 (in), _2 (-y), _3 (+y), _4 (out).
function build_tank(lay, res, nz, zbase)
  (; L, W, Ld_in, Ld_lat, xp₀, xd_out, nx, ny) = lay
  cartesian = G.CartesianDomain(
    mins = (0.0, -W / 2, -H),
    maxs = (L, W / 2, 0.0),
    parts = (nx, ny, nz),
    map = x -> VectorValue(x[1], x[2], G.f_z(x[3] + H, H, nz, zbase) - H),
  )
  G.TankDomain(cartesian;
    structure_domains = [
      G.StructureDomain(L = Lp, W = Bp, x₀ = [xp₀, -Bp / 2, 0.0], domain_symbol = :Γ_plate),
    ],
    damping_zones = [
      G.DampingZone(L = Ld_in, x₀ = [0.0, -W / 2, 0.0], domain_symbol = :Γd_in),
      G.DampingZone(L = xd_out - Ld_in, W = Ld_lat, x₀ = [Ld_in, -W / 2, 0.0], domain_symbol = :Γd_ym),
      G.DampingZone(L = xd_out - Ld_in, W = Ld_lat, x₀ = [Ld_in, W / 2 - Ld_lat, 0.0], domain_symbol = :Γd_yp),
      G.DampingZone(L = lay.Ld_out, x₀ = [xd_out, -W / 2, 0.0], domain_symbol = :Γd_out),
    ],
    resonator_domains = [G.ResonatorDomain(location = loc, trian_symbol = :Γη, delta_symbol = sym)
                         for (loc, sym) in zip(res.locs, res.δsyms)],
  )
end

function damping_coefficients(lay, ws)
  (; W, Ld_in, Ld_out, Ld_lat, xd_out) = lay
  yd = W / 2 - Ld_lat
  s(x) = max(clamp((Ld_in - x[1]) / Ld_in, 0.0, 1.0),
             clamp((x[1] - xd_out) / Ld_out, 0.0, 1.0),
             clamp((abs(x[2]) - yd) / Ld_lat, 0.0, 1.0))
  μ₁(x) = μ₀ * (1.0 - cos(π / 2 * s(x)))
  μ₂(x) = μ₁(x) * ws.kp
  χgen(x) = x[1] < xd_out ? 1.0 : 0.0
  (; μ₁, μ₂, χgen)
end

function _wave(ws, p)
  Dθ = WS.AngularSpreading.DiscreteAngularSpreading(:normal, θ̄, σθ, θmin, θmax, p.nθ)
  Cspec = WS.ContinuousSpectrums.JONSWAP(Hs, Tp)
  Sf = WS.SpectralSpreading.DiscreteSpectralSpreading(
    Cspec, WS.SpectralSampling.UniformSampling(), ws.fmin, ws.fmax, p.nf;
    domain = WS.SpectralSampling.Energy, mess = false,
  )
  s₀ = WS.AiryWaves.AiryState(Sf, Dθ, H)
  sea_state = WS.AiryWaves.AiryState(s₀.spectrum, s₀.spread, s₀.nω, s₀.nθ,
                                     s₀.ω, s₀.k, s₀.θ, s₀.h, p.seed)
  r = WS.realize(sea_state)
  ηin(x, t) = WS.evaluate_η(r, x[1], x[2], t)
  vzin(x, t) = WS.evaluate_w(r, x[1], x[2], x[3], t)
  vin(x, t) = VectorValue(WS.evaluate_u(r, x[1], x[2], x[3], t),
                          WS.evaluate_v(r, x[1], x[2], x[3], t),
                          vzin(x, t))
  (; ηin, vin, vzin)
end

_ramp(t; Tr = 3 * Tp) = t < Tr ? 0.5 * (1 - cos(π * t / Tr)) : 1.0

# ─────────────────────────────────────────────────────────────────────────────
# Solver
# ─────────────────────────────────────────────────────────────────────────────

function main(p; outdir = joinpath(@__DIR__, "..", "data", "tmp_3D_resonator_hl"))
  ws  = wave_scales(p.fmax_fac)
  lay = tank_layout(ws, p.Δ)
  res = resonator_array(lay)
  dz  = damping_coefficients(lay, ws)
  _print_setup(ws, lay, res, p)

  tank = build_tank(lay, res, p.nz, p.zbase)
  wave = _wave(ws, p)

  vin(t)  = x -> _ramp(t) * wave.vin(x, t)
  ηin(t)  = x -> _ramp(t) * wave.ηin(x, t)
  vzin(t) = x -> _ramp(t) * wave.vzin(x, t)

  # DampingZoneBC scales every μ₁ term by βₕ; tmp_3D.jl does not, so divide.
  μ₁ = x -> dz.μ₁(x) / βₕ
  generation_zone(i) = P.DampingZoneBC(
    domain = Symbol("dΓd_$i"), μ₁ = μ₁, μ₂ = dz.μ₂, βₕ = βₕ,
    η_in = ηin, vz_in = vzin,
  )

  fe = PH.FESpaceConfig(order = order, vector_type = Vector{Float64})
  flow = P.PotentialFlow(
    ρw = ρ, g = g, dim = 3, fe = fe, space_domain_symbol = :Ω,
    boundary_conditions = [
      P.PrescribedInletPotentialBC(domain = :dΓin, forcing = vin, quantity = :velocity),
      P.PrescribedInletPotentialBC(domain = :dΓlateral, forcing = vin,
                                   quantity = :velocity, mask = dz.χgen),
      generation_zone(1), generation_zone(2), generation_zone(3),
      P.DampingZoneBC(domain = :dΓd_4, μ₁ = μ₁, μ₂ = dz.μ₂, βₕ = βₕ),
    ],
  )
  free_surface = P.FreeSurface(ρw = ρ, g = g, βₕ = βₕ, dim = 3, fe = fe,
                               space_domain_symbol = :Γκ)
  plate = P.KirchhoffLovePlate(
    E = E, ν = ν, hb = hb, ρ = ρ, ρb = ρb, g = g,
    fe = PH.FESpaceConfig(order = order, vector_type = Vector{Float64},
                          γ = 1.0 * order * (order + 1)),
    space_domain_symbol = :Γη,
  )
  resonators = P.resonator_array(res.N, res.M, res.K, res.C; ρw = ρ,
                                 delta_domain_symbols = res.δsyms,
                                 fe = PH.FESpaceConfig(order = order, vector_type = Vector{Float64}))

  config  = S.TimeDomainConfig(t₀ = 0.0, tf = p.T)
  tconfig = S.TimeConfig(Δt = p.Δt, t₀ = 0.0, tf = p.T, scheme = :newmark,
                         γ = γₜ, β = βₜ)
  problem = S.build_problem(tank, P.PhysicsParameters[flow, free_surface, plate, resonators],
                            config; tconfig = tconfig)
  trians = S.get_triangulations(problem)
  dom    = S.get_integration_domains(problem)
  Γκ, Γη = trians[:Γκ], trians[:Γη]
  δ = [dom[s] for s in res.δsyms]
  result = S.simulate(problem, tconfig)

  mkpath(outdir)
  open(joinpath(outdir, "resonators.csv"), "w") do f
    println(f, "x,y,M,K,C,omega_r")
    for i in 1:res.N
      println(f, res.locs[i][1], ",", res.locs[i][2], ",", res.M[i], ",", res.K[i], ",",
              res.C[i], ",", res.ωr[i])
    end
  end
  io = open(joinpath(outdir, "resonator_timeseries.txt"), "w")
  println(io, "# t  " * join(["eta_$i q_$i" for i in 1:res.N], "  "))
  nout = max(1, round(Int, p.Δt_out / p.Δt))
  if p.vtk_output
    writevtk(Γη, joinpath(outdir, "plate_triangulation"))
    pvdη = createpvd(joinpath(outdir, "solution_η"))
    pvdκ = createpvd(joinpath(outdir, "solution_κ"))
  end
  for (t, x) in result.solution
    istep = round(Int, t / p.Δt)
    κₕ, ηₕ = x[2], x[3]
    ηp = [sum(δ[i](ηₕ)) for i in 1:res.N]
    qp = [get_free_dof_values(x[3 + i])[1] for i in 1:res.N]
    println(io, t, "  ", join(("$(ηp[i]) $(qp[i])" for i in 1:res.N), "  "))
    println("t = $t: max|η(pᵢ)| = $(maximum(abs, ηp)), max|qᵢ| = $(maximum(abs, qp))")
    if p.vtk_output && istep % nout == 0
      # N-component constant fields; component i follows the row order of resonators.csv
      q_all = CellField(VectorValue(qp...), Γη)
      η_res = CellField(VectorValue(ηp...), Γη)
      pvdη[t] = createvtk(Γη, joinpath(outdir, "solution_η_$istep"),
                          cellfields = ["η" => ηₕ, "q" => q_all, "eta_r" => η_res]; order = order)
      pvdκ[t] = createvtk(Γκ, joinpath(outdir, "solution_κ_$istep"),
                          cellfields = ["κ" => κₕ]; order = order)
    end
  end
  close(io)
  if p.vtk_output
    savepvd(pvdη)
    savepvd(pvdκ)
  end

  return nothing
end

function _print_setup(ws, lay, res, p)
  println("── Tank ───────────────────────────────────────────")
  println("  L × W × H = $(lay.L) × $(lay.W) × $H m, nx × ny × nz = $(lay.nx) × $(lay.ny) × $(p.nz), Δ = $(p.Δ) m")
  println("  plate x ∈ [$(lay.xp₀), $(lay.xp₀ + Lp)], y ∈ [$(-Bp/2), $(Bp/2)]")
  println("── Resonators ─────────────────────────────────────")
  println("  $(nrx) × $(nry) = $(res.N) resonators, M = $(round(res.M[1], digits=1)) kg each, ζr = $ζr")
  println("  ωr graded in x: $(round(ωr_lo, digits=3)) → $(round(ωr_hi, digits=3)) rad/s ",
          "(ωp = $(round(ωp, digits=3)) rad/s), ",
          "K = $(round(minimum(res.K), digits=1)) → $(round(maximum(res.K), digits=1)) N/m")
end

# Smoke test only (2 steps on a coarse mesh); DEV/PROD are real simulations.
@timeit to "smoke" main(SMOKE)
@timeit to "dev"  main(DEV)
# @timeit to "prod" main(PROD)

show(to)

end
