module TMP3D

# 3D time-domain hydroelastic response of a floating Kirchhoff-Love plate in a
# short-crested JONSWAP sea, built on the `TankDomain{3}` geometry of
# HydroElasticFEM.jl (branch dev/tankdomain3d).
#
# Tank layout (top view, x = main wave direction, free surface at z = 0):
#
#   y = +W/2 ┌──────┬──────────── lateral zone (+y) ────────────┬────────┐
#            │      │                                           │        │
#            │inlet │  G_up     ┌───────┐    G_down             │ outlet │
#            │ zone │ ────────▶ │ plate │ ─────────▶            │  zone  │
#            │      │           └───────┘                       │        │
#            │      │                                           │        │
#   y = -W/2 └──────┴──────────── lateral zone (-y) ────────────┴────────┘
#            0    Ld_in                                        xd_out     L
#
# - Inlet and lateral zones relax (κ, ∂ₙϕ) towards the incident sea: they
#   generate the incident field and absorb scattered / radiated waves.
# - The outlet zone relaxes towards zero (pure absorber).
# - The incident normal flux is prescribed on the inlet wall and on the
#   lateral walls upstream of the outlet zone (the spread sea is not
#   compatible with a no-flux lateral wall).

import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import WaveSpec as WS
using Gridap
using Gridap.ODEs
using TimerOutputs

to = TimerOutput()

# ─────────────────────────────────────────────────────────────────────────────
# Physical parameters
# ─────────────────────────────────────────────────────────────────────────────
H  = 30.0      # Water depth [m]
ρ  = 1000.0    # Water density [kg/m³]
g  = 9.81      # Gravity [m/s²]

# Plate (Yago-type VLFS material, 60 × 50 m)
Lp = 60.0      # Plate length along x [m]
Bp = 50.0      # Plate width along y [m]
E  = 11.9e9
ρb = 256.25
ν  = 0.13
hb = 2.0
D  = E * hb^3 / 12 / (1 - ν^2)
d₀ = ρb * hb / ρ
C  = P.build_KL_tensor(E, ν, hb, ρ)

# Sea state
Hs = 3.0       # Significant wave height [m]
Tp = 10.0      # Peak period [s]
θ̄  = 0.0       # Mean direction [rad]
σθ = π / 6     # Directional standard deviation [rad]
θmin, θmax = -π / 4, π / 4

# Spectral band kept in the simulation, as multiples of fp.  For JONSWAP the
# energy below 0.7 fp is < 0.5 % and above 2.5 fp ≈ 2 % of m₀.
fmin_fac = 0.7

# Damping law (Colomés et al. 2023): μ₁ = μ₀ (1 - cos(π/2 s)), μ₂ = μ₁ kₚ,
# with s ∈ [0, 1] the normalised depth into the zone.
μ₀ = 2.5

# ─────────────────────────────────────────────────────────────────────────────
# Discretisation presets
# ─────────────────────────────────────────────────────────────────────────────
#   Δ        horizontal cell size (Δx = Δy); Q2 elements need ≥ 5 cells per λmin
#   nz,zbase vertical cells and grading base of G.f_z (fine near the surface)
#   fmax_fac upper spectral cut-off (× fp), sets λmin = λ(fmax)
#   nf, nθ   frequency / direction edges (bins = n - 1)
#   Δt       Newmark step: ≤ Tmin/25 (dev) … Tmin/40 (production)
#   T        simulated time: ~20 Tp transient + statistics window
#   Δt_out   VTK output interval
const WARMUP = (; Δ = 10.0, nz = 1, zbase = 2.0, fmax_fac = 2.0, nf = 40, nθ = 11,
               Δt = 0.2, T = 0.2, Δt_out = 0.2, seed = 1234, vtk_output = true)
const DEV = (; Δ = 10.0, nz = 4, zbase = 2.0, fmax_fac = 2.0, nf = 40, nθ = 11,
               Δt = 0.2, T = 4.0, Δt_out = 1.0, seed = 1234, vtk_output = true)
const PROD = (; Δ = 5.0, nz = 5, zbase = 2.0, fmax_fac = 2.5, nf = 50, nθ = 11,
                Δt = 0.1, T = 1400.0, Δt_out = 0.5, seed = 1234, vtk_output = true)

# ─────────────────────────────────────────────────────────────────────────────
# Wave scales and tank layout
# ─────────────────────────────────────────────────────────────────────────────

"""
    wave_scales(fmax_fac)

Peak, longest and shortest wavelengths of the simulated band at depth `H`.
"""
function wave_scales(fmax_fac)
  fp = 1 / Tp
  fmin, fmax = fmin_fac * fp, fmax_fac * fp
  k(f) = WS.AiryWaves.solve_wavenumber(2π * f, H)
  kp = k(fp)
  (; fp, fmin, fmax, kp, λp = 2π / kp, λmax = 2π / k(fmin), λmin = 2π / k(fmax),
     Tmin = 1 / fmax)
end

# Round a length up to a whole number of cells, so every zone edge falls on a
# cell face and the centroid masks of TankDomain select whole cells.
_snap(x, Δ) = Δ * ceil(x / Δ - 1e-9)

"""
    tank_layout(ws, Δ)

Zone lengths derived from the sea state (Tp = 10 s, H = 30 m gives
λp ≈ 137 m, λmax ≈ 221 m):

- `Ld_in  ≈ λmax`   inlet generation/absorption zone
- `Ld_out ≈ 2 λp`   outlet absorber (≥ 1.2 λmax)
- `Ld_lat ≈ λp`     lateral generation/absorption zones
- `G_up = G_down ≈ λp`  free water between zones and plate
- `G_lat ≈ λp/2`    free water between plate side and lateral zone
"""
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

"""
    build_tank(lay, nz, zbase)

3D `TankDomain` with free surface at z = 0, graded vertical mesh, the plate as
a `StructureDomain` and four non-overlapping damping zones.
"""
function build_tank(lay, nz, zbase)
  (; L, W, Ld_in, Ld_lat, xp₀, xd_out, nx, ny) = lay
  cartesian = G.CartesianDomain(
    mins = (0.0, -W / 2, -H),
    maxs = (L, W / 2, 0.0),
    parts = (nx, ny, nz),
    # G.f_z grades [0, H] towards z = H; shift so the surface stays at z = 0
    # (WaveSpec convention: z ∈ [-h, 0]).
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

"""
    damping_coefficients(lay, ws)

Smooth relaxation coefficients `μ₁`, `μ₂` (zero outside the zones) and the
indicator `χgen` of the region where the target is the incident field.
"""
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

function _wave(ws, p, Ω, old_approach::Bool)
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
  # Fix the random phases so DEV and PROD runs see the same realisation
  s₀ = WS.AiryWaves.AiryState(Sf, Dθ, H)
  sea_state = WS.AiryWaves.AiryState(s₀.spectrum, s₀.spread, s₀.nω, s₀.nθ,
                                     s₀.ω, s₀.k, s₀.θ, s₀.h, p.seed)
  if old_approach
    coords = unique(reduce(vcat, get_cell_coordinates(Ω)))
    x_coords = unique([c[1] for c in coords])
    y_coords = unique([c[2] for c in coords])
    z_coords = unique([c[3] for c in coords])
    t_vals = 0:p.Δt:p.T
    println("Generating interpolable sea...")
    sea = WS.AiryWaves.generate_interpolable_sea(sea_state, x_coords, y_coords, z_coords, t_vals)
    println("Interpolable sea generated.")
    ηin_old(x, t) = sea[:η](x[1], x[2], x[3], t)
    ϕin_old(x, t) = sea[:ϕ](x[1], x[2], x[3], t)
    vzin_old(x, t) = sea[:w](x[1], x[2], x[3], t)
    vin_old(x, t) = VectorValue(sea[:u](x[1], x[2], x[3], t), sea[:v](x[1], x[2], x[3], t), vzin_old(x, t))
    return (; ηin = ηin_old, ϕin = ϕin_old, vin = vin_old, vzin = vzin_old)
  else
    r = WS.realize(sea_state)
    ηin_new(x, t) = WS.evaluate_η(r, x[1], x[2], t)
    ϕin_new(x, t) = WS.evaluate_ϕ(r, x[1], x[2], x[3], t)
    vzin_new(x, t) = WS.evaluate_w(r, x[1], x[2], x[3], t)
    vin_new(x, t) = VectorValue(WS.evaluate_u(r, x[1], x[2], x[3], t),
                                WS.evaluate_v(r, x[1], x[2], x[3], t),
                                vzin_new(x, t))
    return (; ηin = ηin_new, ϕin = ϕin_new, vin = vin_new, vzin = vzin_new)
  end
end

# Smooth start of the forcing over 3 Tp to avoid an impulsive transient
_ramp(t; Tr = 3 * Tp) = t < Tr ? 0.5 * (1 - cos(π * t / Tr)) : 1.0

# ─────────────────────────────────────────────────────────────────────────────
# Solver
# ─────────────────────────────────────────────────────────────────────────────

function main(p; old_approach::Bool = false, outdir = joinpath(@__DIR__, "..", "data", "tmp_3D"))
  ws  = wave_scales(p.fmax_fac)
  lay = tank_layout(ws, p.Δ)
  dz  = damping_coefficients(lay, ws)
  _print_setup(ws, lay, p)

  # Geometry: every triangulation and measure comes from the TankDomain
  #   Ω  fluid volume      Γκ free surface incl. damping zones (κ field)
  #   Γη plate (η field)   Λη plate interior skeleton (C/DG terms)
  order = 2
  tank   = build_tank(lay, p.nz, p.zbase)
  model  = G.build_model(tank)
  trians = G.build_triangulations(tank, model)
  dom    = G.get_integration_domains(trians; degree = 2 * order)
  Ω, Γκ, Γη = trians[:Ω], trians[:Γκ], trians[:Γη]
  dΩ, dΓκ, dΓη, dΓin = dom[:dΩ], dom[:dΓκ], dom[:dΓη], dom[:dΓin]
  dΛη, nΛη = dom[:dΛη], dom[:n_Λ_η]
  nΓin = get_normal_vector(trians[:Γin])
  Γlat = trians[:Γlateral]
  dΓlat, nΓlat = Measure(Γlat, 2 * order), get_normal_vector(Γlat)
  dΓgen = [Measure(trians[s], 2 * order) for s in (:Γd_in, :Γd_ym, :Γd_yp)]
  if p.vtk_output
    mkpath(outdir)
    writevtk(Γη, joinpath(outdir, "plate_triangulation"))
    writevtk(Γκ, joinpath(outdir, "surface_triangulation"),
             cellfields = ["μ₁" => CellField(dz.μ₁, Γκ)])
  end

  wave = _wave(ws, p, Ω, old_approach)
  vin(t)  = (r = _ramp(t); x -> r * wave.vin(x, t))
  ηd(t)   = (r = _ramp(t); x -> r * dz.μ₂(x) * wave.ηin(x, t))
  ∇ₙϕd(t) = (r = _ramp(t); x -> r * dz.μ₁(x) * wave.vzin(x, t))

  # Spaces
  reffe = ReferenceFE(lagrangian, Float64, order)
  Vϕ = TestFESpace(Ω, reffe; conformity = :H1)
  Vκ = TestFESpace(Γκ, reffe; conformity = :H1)
  Vη = TestFESpace(Γη, reffe; conformity = :H1)
  Y = MultiFieldFESpace([Vϕ, Vκ, Vη])
  X = MultiFieldFESpace([TrialFESpace(Vϕ), TrialFESpace(Vκ), TrialFESpace(Vη)])

  # Coefficients
  βₕ = 0.5
  γₜ = 0.5
  βₜ = 0.25
  αₕ = γₜ / (βₜ * p.Δt) / g * (1.0 - βₕ) / βₕ
  γ  = 1.0 * order * (order + 1) / dom[:h_η]   # C/DG penalty, h = Δ

  # Weak form
  ∇ₙ(ϕ) = ∇(ϕ) ⋅ VectorValue(0.0, 0.0, 1.0)
  m(t, (ϕₜₜ, κₜₜ, ηₜₜ), (w, u, v)) = ∫( d₀ * ηₜₜ * v )dΓη
  c(t, (ϕₜ, κₜ, ηₜ), (w, u, v)) =
    ∫( βₕ * ϕₜ * (u + αₕ * w) - κₜ * w )dΓκ +
    ∫( ϕₜ * v - ηₜ * w )dΓη
  k(t, (ϕ, κ, η), (w, u, v)) =
    ∫( ∇(w) ⋅ ∇(ϕ) )dΩ +
    # free surface + relaxation (μ₁ = μ₂ = 0 outside the damping zones)
    ∫( βₕ * (u + αₕ * w) * (g * κ) - dz.μ₂ * κ * w + dz.μ₁ * ∇ₙ(ϕ) * (u + αₕ * w) )dΓκ +
    # plate: hydrostatics + Kirchhoff-Love bending
    ∫( v * (g * η) + ∇∇(v) ⊙ (C ⊙ ∇∇(η)) )dΓη +
    # C/DG continuity of the plate slope across interior facets
    ∫( -jump(∇(v)) ⊙ (mean(C ⊙ ∇∇(η)) ⋅ nΛη.⁺) -
       (mean(C ⊙ ∇∇(v)) ⋅ nΛη.⁺) ⊙ jump(∇(η)) +
       D / ρ * γ * jump(∇(v)) ⊙ jump(∇(η)) )dΛη
  b(t, (w, u, v)) =
    ∫( w * (vin(t) ⋅ nΓin) )dΓin +
    ∫( w * dz.χgen * (vin(t) ⋅ nΓlat) )dΓlat +
    sum(∫( ∇ₙϕd(t) * (u + αₕ * w) - ηd(t) * w )dΓ for dΓ in dΓgen)

  op = TransientLinearFEOperator((k, c, m), b, X, Y; constant_forms = (true, true, true))
  ode_solver = Newmark(LUSolver(), p.Δt, γₜ, βₜ)

  # Quiescent initial state (forcing is ramped)
  x₀ = interpolate_everywhere([0.0, 0.0, 0.0], X(0.0))
  v₀ = interpolate_everywhere([0.0, 0.0, 0.0], X(0.0))
  a₀ = interpolate_everywhere([0.0, 0.0, 0.0], X(0.0))

  xₜ = solve(ode_solver, op, 0.0, p.T, (x₀, v₀, a₀))

  nout = max(1, round(Int, p.Δt_out / p.Δt))
  if p.vtk_output
    pvdη = createpvd(joinpath(outdir, "solution_η"))
    pvdκ = createpvd(joinpath(outdir, "solution_κ"))
    pvdϕ = createpvd(joinpath(outdir, "solution_ϕ"))
  end
  for (t, x) in xₜ
    istep = round(Int, t / p.Δt)
    println("solution at t = ", t)
    if p.vtk_output && istep % nout == 0
      ϕₕ, κₕ, ηₕ = x
      pvdη[t] = createvtk(Γη, joinpath(outdir, "solution_η_$istep"), cellfields = ["η" => ηₕ]; order = order)
      pvdκ[t] = createvtk(Γκ, joinpath(outdir, "solution_κ_$istep"),
                          cellfields = ["κ" => κₕ, "ηin" => CellField(x -> _ramp(t) * wave.ηin(x, t), Γκ)]; order = order)
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
@timeit to "warm_up" main(WARMUP)
@timeit to "dev"     main(DEV)
# @timeit to "prod"    main(PROD)

show(to)

end
