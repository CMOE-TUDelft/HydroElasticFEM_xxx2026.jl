module TMP3DHL

# Same problem as tmp_3D.jl (floating Kirchhoff-Love plate in a short-crested
# JONSWAP sea, 3D time domain), but assembled with the high-level API of
# HydroElasticFEM.jl: build_problem + simulate, PotentialFlow boundary
# conditions for the inlet/lateral flux and the relaxation zones.
#
# Tank layout (top view, x = main wave direction, free surface at z = 0):
#
#   y = +W/2 ┌──────┬──────────── lateral zone (+y) ────────────┬────────┐
#            │inlet │  G_up     ┌───────┐    G_down             │ outlet │
#            │ zone │ ────────▶ │ plate │ ─────────▶            │  zone  │
#   y = -W/2 └──────┴──────────── lateral zone (-y) ────────────┴────────┘
#            0    Ld_in                                        xd_out     L

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
H  = 30.0      # Water depth [m]
ρ  = 1000.0    # Water density [kg/m³]
g  = 9.81      # Gravity [m/s²]

# Plate (Yago-type VLFS material, 60 × 50 m)
Lp = 60.0
Bp = 50.0
E  = 11.9e9
ρb = 256.25
ν  = 0.13
hb = 2.0

# Sea state
Hs = 3.0
Tp = 10.0
θ̄  = 0.0
σθ = π / 6
θmin, θmax = -π / 4, π / 4
fmin_fac = 0.7

# Damping law: μ₁ = μ₀ (1 - cos(π/2 s)), μ₂ = μ₁ kₚ, s ∈ [0, 1] depth into the zone.
μ₀ = 2.5

# Free-surface stabilisation and Newmark parameters
βₕ = 0.5
γₜ = 0.5
βₜ = 0.25
order = 2

# ─────────────────────────────────────────────────────────────────────────────
# Discretisation presets (see tmp_3D.jl for the meaning of each entry)
# ─────────────────────────────────────────────────────────────────────────────
const WARM = (; Δ = 10.0, nz = 1, zbase = 2.0, fmax_fac = 2.0, nf = 40, nθ = 11,
               Δt = 0.2, T = 0.2, Δt_out = 0.2, seed = 1234, vtk_output = true)
const DEV = (; Δ = 10.0, nz = 4, zbase = 2.0, fmax_fac = 2.0, nf = 40, nθ = 11,
               Δt = 0.2, T = 4.0, Δt_out = 1.0, seed = 1234, vtk_output = true)
const PROD = (; Δ = 5.0, nz = 5, zbase = 2.0, fmax_fac = 2.5, nf = 50, nθ = 11,
                Δt = 0.1, T = 1400.0, Δt_out = 0.5, seed = 1234, vtk_output = true)

# ─────────────────────────────────────────────────────────────────────────────
# Wave scales and tank layout
# ─────────────────────────────────────────────────────────────────────────────

function wave_scales(fmax_fac)
  fp = 1 / Tp
  fmin, fmax = fmin_fac * fp, fmax_fac * fp
  k(f) = WS.AiryWaves.solve_wavenumber(2π * f, H)
  kp = k(fp)
  (; fp, fmin, fmax, kp, λp = 2π / kp, λmax = 2π / k(fmin), λmin = 2π / k(fmax),
     Tmin = 1 / fmax)
end

# Round up to a whole number of cells so zone edges fall on cell faces.
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

# Damping zones are numbered in this order: :dΓd_1 (in), _2 (-y), _3 (+y), _4 (out).
function build_tank(lay, nz, zbase)
  (; L, W, Ld_in, Ld_lat, xp₀, xd_out, nx, ny) = lay
  cartesian = G.CartesianDomain(
    mins = (0.0, -W / 2, -H),
    maxs = (L, W / 2, 0.0),
    parts = (nx, ny, nz),
    # G.f_z grades [0, H] towards z = H; shift so the surface stays at z = 0.
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

# ─────────────────────────────────────────────────────────────────────────────
# Incident sea
# ─────────────────────────────────────────────────────────────────────────────

function _wave(ws, p)
  Dθ = WS.AngularSpreading.DiscreteAngularSpreading(:normal, θ̄, σθ, θmin, θmax, p.nθ)
  Cspec = WS.ContinuousSpectrums.JONSWAP(Hs, Tp)
  Sf = WS.SpectralSpreading.DiscreteSpectralSpreading(
    Cspec,
    WS.SpectralSampling.UniformSampling(),
    ws.fmin,
    ws.fmax,
    p.nf;
    domain = WS.SpectralSampling.Energy, mess = false,
  )
  # Fixed random phases so DEV and PROD runs see the same realisation
  s₀ = WS.AiryWaves.AiryState(Sf, Dθ, H)
  sea_state = WS.AiryWaves.AiryState(s₀.spectrum, s₀.spread, s₀.nω, s₀.nθ,
                                     s₀.ω, s₀.k, s₀.θ, s₀.h, p.seed)
  r = WS.realize(sea_state)
  ηin(x, t) = WS.evaluate_η(r, x[1], x[2], t)
  vzin(x, t) = WS.evaluate_w(r, x[1], x[2], x[3], t)
  vin(x, t) = VectorValue(WS.evaluate_u(r, x[1], x[2], x[3], t),
                          WS.evaluate_v(r, x[1], x[2], x[3], t),
                          vzin(x, t))
  (; ηin, vzin, vin)
end

# Smooth start of the forcing over 3 Tp to avoid an impulsive transient
_ramp(t; Tr = 3 * Tp) = t < Tr ? 0.5 * (1 - cos(π * t / Tr)) : 1.0

# ─────────────────────────────────────────────────────────────────────────────
# Solver
# ─────────────────────────────────────────────────────────────────────────────

function main(p; outdir = joinpath(@__DIR__, "..", "data", "tmp_3D_hl"))
  ws  = wave_scales(p.fmax_fac)
  lay = tank_layout(ws, p.Δ)
  dz  = damping_coefficients(lay, ws)
  _print_setup(ws, lay, p)

  tank = build_tank(lay, p.nz, p.zbase)
  wave = _wave(ws, p)

  # Time-dependent forcing, resolved as t -> x -> value
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
  # C/DG penalty order*(order+1)/h with h = Δ
  plate = P.KirchhoffLovePlate(
    E = E, ν = ν, hb = hb, ρ = ρ, ρb = ρb, g = g,
    fe = PH.FESpaceConfig(order = order, vector_type = Vector{Float64},
                          γ = 1.0 * order * (order + 1)),
    space_domain_symbol = :Γη,
  )

  config  = S.TimeDomainConfig(t₀ = 0.0, tf = p.T)
  tconfig = S.TimeConfig(Δt = p.Δt, t₀ = 0.0, tf = p.T, scheme = :newmark,
                         γ = γₜ, β = βₜ)
  problem = S.build_problem(tank, P.PhysicsParameters[flow, free_surface, plate],
                            config; tconfig = tconfig)
  trians = S.get_triangulations(problem)
  Ω, Γκ, Γη = trians[:Ω], trians[:Γκ], trians[:Γη]

  if p.vtk_output
    mkpath(outdir)
    writevtk(Γη, joinpath(outdir, "plate_triangulation"))
    writevtk(Γκ, joinpath(outdir, "surface_triangulation"),
             cellfields = ["μ₁" => CellField(dz.μ₁, Γκ)])
  end

  result = S.simulate(problem, tconfig)

  nout = max(1, round(Int, p.Δt_out / p.Δt))
  if p.vtk_output
    pvdη = createpvd(joinpath(outdir, "solution_η"))
    pvdκ = createpvd(joinpath(outdir, "solution_κ"))
    pvdϕ = createpvd(joinpath(outdir, "solution_ϕ"))
  end
  for (t, x) in result.solution
    istep = round(Int, t / p.Δt)
    println("solution at t = ", t)
    if p.vtk_output && istep % nout == 0
      ϕₕ, κₕ, ηₕ = x
      pvdη[t] = createvtk(Γη, joinpath(outdir, "solution_η_$istep"), cellfields = ["η" => ηₕ]; order = order)
      pvdκ[t] = createvtk(Γκ, joinpath(outdir, "solution_κ_$istep"),
                          cellfields = ["κ" => κₕ, "ηin" => CellField(ηin(t), Γκ)]; order = order)
      pvdϕ[t] = createvtk(Ω, joinpath(outdir, "solution_ϕ_$istep"), cellfields = ["ϕ" => ϕₕ])
    end
  end
  if p.vtk_output
    savepvd(pvdη)
    savepvd(pvdκ)
    savepvd(pvdϕ)
  end

  return nothing
end

function _print_setup(ws, lay, p)
  println("── Sea state ──────────────────────────────────────")
  println("  λp = $(round(ws.λp, digits=1)) m, λmax = $(round(ws.λmax, digits=1)) m (f = $(ws.fmin) Hz), ",
          "λmin = $(round(ws.λmin, digits=1)) m (f = $(ws.fmax) Hz)")
  println("── Tank layout ────────────────────────────────────")
  println("  L × W × H = $(lay.L) × $(lay.W) × $H m;  Ld_in = $(lay.Ld_in), Ld_out = $(lay.Ld_out), ",
          "Ld_lat = $(lay.Ld_lat), G_up = $(lay.G_up), G_down = $(lay.G_down), G_lat = $(lay.G_lat)")
  println("  plate x ∈ [$(lay.xp₀), $(lay.xp₀ + Lp)], y ∈ [$(-Bp/2), $(Bp/2)]")
  println("── Discretisation ─────────────────────────────────")
  println("  nx × ny × nz = $(lay.nx) × $(lay.ny) × $(p.nz),  Δ = $(p.Δ) m ",
          "($(round(ws.λmin / p.Δ, digits=1)) Q2 cells per λmin),  Δt = $(p.Δt) s ",
          "(Tmin/Δt = $(round(ws.Tmin / p.Δt, digits=1)))")
  ws.λmin / p.Δ < 5 && @warn "Fewer than 5 Q2 cells per shortest wavelength."
end

# Warm-up (compilation), development run, production run
@timeit to "warm_up" main(WARM)
@timeit to "dev"     main(DEV)
# @timeit to "prod"    main(PROD)

show(to)

end
