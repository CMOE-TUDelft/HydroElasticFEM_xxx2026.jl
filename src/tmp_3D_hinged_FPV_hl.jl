module TMP3DHingedFPVHL

# Same problem as tmp_3D_hinged_FPV.jl (14 × 14 hinged floaters, Nymphaea
# Aurora-type FPV array, short-crested JONSWAP sea, 3D time domain), assembled
# with the high-level API of HydroElasticFEM.jl: `G.hinge_grid` for the line
# hinges, `KirchhoffLovePlate`, build_problem + simulate.  See that file for
# the site data and modelling assumptions.
#
# Top view (x = main wave direction, free surface at z = 0):
#
#   y = +W/2 ┌──────┬──────────── lateral zone (+y) ────────────┬────────┐
#            │inlet │  G_up     ┌─┬─┬─┬─┬─┐        G_down       │ outlet │
#            │ zone │ ────────▶ ├─┼─┼─┼─┼─┤ ─────────▶          │  zone  │
#            │      │           └─┴─┴─┴─┴─┘  hinged floaters    │        │
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
# Site and floater array
# ─────────────────────────────────────────────────────────────────────────────
H  = 25.0
ρ  = 1025.0
g  = 9.81

a   = 7.0      # Floater side length [m]
nfx = 14
nfy = 14
Lp  = nfx * a
Bp  = nfy * a

# Equivalent floater: areal mass ms and bending stiffness D (see tmp_3D_hinged_FPV.jl)
ms  = 100.0
D   = 1.0e8
E   = 70.0e9
ν   = 0.3
hb  = cbrt(12 * (1 - ν^2) * D / E)
ρb  = ms / hb                            # so that ρb hb / ρ = ms / ρ
kᵣ  = 0.0      # Hinge rotational stiffness per unit length [N·m/rad/m]; 0 = free hinge

# Sea state
Hs = 2.5
Tp = 7.5
γJ = 3.3
θ̄  = 0.0
σθ = deg2rad(25.0)
θmin, θmax = -π / 3, π / 3
fmin_fac = 0.7

μ₀ = 2.5

βₕ = 0.5
γₜ = 0.5
βₜ = 0.25
order = 2

# ─────────────────────────────────────────────────────────────────────────────
# Discretisation presets (see tmp_3D_hinged_FPV.jl for the meaning of each entry)
# ─────────────────────────────────────────────────────────────────────────────
const WARM = (; Δc = 14.0, Δf = 7.0, nz = 1, zbase = 2.0, fmax_fac = 1.2, nf = 20, nθ = 7,
               Δt = 0.3, T = 0.6, Δt_out = 0.3, seed = 1234, vtk_output = true)
const DEV = (; Δc = 7.0, Δf = 3.5, nz = 4, zbase = 2.0, fmax_fac = 1.5, nf = 40, nθ = 11,
               Δt = 0.15, T = 3.0, Δt_out = 0.6, seed = 1234, vtk_output = true)
const PROD = (; Δc = 3.5, Δf = 1.75, nz = 5, zbase = 2.0, fmax_fac = 2.2, nf = 60, nθ = 15,
                Δt = 0.08, T = 1200.0, Δt_out = 0.4, seed = 1234, vtk_output = true)

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

_snap(x, Δ) = Δ * ceil(x / Δ - 1e-9)
_is_multiple(x, Δ) = abs(x / Δ - round(x / Δ)) < 1e-9

# Face coordinates of a 1D grid from `segments = [(length, Δ), …]`
function _faces(x₀, segments)
  xf = [x₀]
  for (len, Δ) in segments
    _is_multiple(len, Δ) || error("Segment length $len is not a multiple of Δ = $Δ.")
    n = round(Int, len / Δ)
    append!(xf, xf[end] .+ Δ .* (1:n))
  end
  xf
end

# Piecewise-linear map from a uniform grid on [xf[1], xf[end]] to the faces `xf`
function _face_map(xf)
  n = length(xf) - 1
  x₀, span = xf[1], xf[end] - xf[1]
  ξ -> begin
    s = (ξ - x₀) / span * n
    i = clamp(floor(Int, s), 0, n - 1)
    xf[i + 1] + (s - i) * (xf[i + 2] - xf[i + 1])
  end
end

# The fine-mesh block covers the array plus one floater `a` of margin.
function tank_layout(ws, p)
  (; Δc, Δf) = p
  for (name, x) in (("a", a), ("Lp", Lp), ("Bp", Bp))
    _is_multiple(x, Δf) || error("$name = $x m must be a multiple of Δf = $(Δf) m.")
  end
  mf = a
  Ld_in  = _snap(ws.λmax, Δc)
  Ld_out = _snap(2 * ws.λp, Δc)
  Ld_lat = _snap(ws.λp, Δc)
  G_up   = mf + _snap(ws.λp - mf, Δc)
  G_down = mf + _snap(ws.λp - mf, Δc)
  G_lat  = mf + _snap(ws.λp / 2 - mf, Δc)

  L  = Ld_in + G_up + Lp + G_down + Ld_out
  W  = Bp + 2 * (G_lat + Ld_lat)
  xp₀ = Ld_in + G_up
  yp₀ = -Bp / 2
  xd_out = L - Ld_out

  xf = _faces(0.0, [(Ld_in, Δc), (G_up - mf, Δc), (Lp + 2mf, Δf), (G_down - mf, Δc), (Ld_out, Δc)])
  yf = _faces(-W / 2, [(Ld_lat, Δc), (G_lat - mf, Δc), (Bp + 2mf, Δf), (G_lat - mf, Δc), (Ld_lat, Δc)])
  (; L, W, Ld_in, Ld_out, Ld_lat, G_up, G_down, G_lat, xp₀, yp₀, xd_out, xf, yf,
     nx = length(xf) - 1, ny = length(yf) - 1)
end

# Damping zones are numbered :dΓd_1 (in), _2 (-y), _3 (+y), _4 (out); the hinge
# lines of the array are one JointLineDomain (:dΛh, :n_Λh).
function build_tank(lay, nz, zbase)
  (; L, W, Ld_in, Ld_lat, xp₀, yp₀, xd_out, nx, ny) = lay
  mx, my = _face_map(lay.xf), _face_map(lay.yf)
  cartesian = G.CartesianDomain(
    mins = (0.0, -W / 2, -H),
    maxs = (L, W / 2, 0.0),
    parts = (nx, ny, nz),
    map = x -> VectorValue(mx(x[1]), my(x[2]), G.f_z(x[3] + H, H, nz, zbase) - H),
  )
  G.TankDomain(cartesian;
    structure_domains = [
      G.StructureDomain(L = Lp, W = Bp, x₀ = [xp₀, yp₀, 0.0], domain_symbol = :Γ_array),
    ],
    joint_domains = [
      G.hinge_grid(x₀ = [xp₀, yp₀, 0.0], a = a, nfx = nfx, nfy = nfy,
                   domain_symbol = :dΛh, normal_symbol = :n_Λh),
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
  Cspec = WS.ContinuousSpectrums.JONSWAP(Hs, Tp, γJ)
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

function main(p; outdir = joinpath(@__DIR__, "..", "data", "tmp_3D_hinged_FPV_hl"))
  ws  = wave_scales(p.fmax_fac)
  lay = tank_layout(ws, p)
  dz  = damping_coefficients(lay, ws)
  _print_setup(ws, lay, p)

  tank = build_tank(lay, p.nz, p.zbase)
  wave = _wave(ws, p)

  vin(t)  = x -> _ramp(t) * wave.vin(x, t)
  ηin(t)  = x -> _ramp(t) * wave.ηin(x, t)
  vzin(t) = x -> _ramp(t) * wave.vzin(x, t)

  # DampingZoneBC scales every μ₁ term by βₕ; tmp_3D_hinged_FPV.jl does not, so divide.
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
  # The hinge facets are split off the C/DG skeleton; a free hinge needs no joint entry.
  joints = kᵣ > 0 ? [P.JointRotationalSpring(:dΛh, :n_Λh, kᵣ / ρ)] : P.JointRotationalSpring[]
  plate = P.KirchhoffLovePlate(
    E = E, ν = ν, hb = hb, ρ = ρ, ρb = ρb, g = g, joints = joints,
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
  dom    = S.get_integration_domains(problem)
  Ω, Γκ, Γη = trians[:Ω], trians[:Γκ], trians[:Γη]
  Λh = trians[:Λ_joints][1]
  dΓη, dΛh = dom[:dΓη], dom[:dΛh]

  Aarr = sum(∫(1.0)dΓη)
  Lh   = sum(∫(1.0)dΛh)
  println("  array area = $(round(Aarr, digits=1)) m², hinge length = $(round(Lh, digits=1)) m")

  # Per-cell floater index (i, j); `checker` alternates between neighbours for renderings
  xcη = map(xs -> sum(xs) / length(xs), collect(get_cell_coordinates(Γη)))
  ifl = [floor(Int, (x[1] - lay.xp₀) / a) for x in xcη]
  jfl = [floor(Int, (x[2] - lay.yp₀) / a) for x in xcη]
  floater = CellField(Float64.(ifl .+ nfx .* jfl), Γη)
  checker = CellField(Float64.(mod.(ifl .+ jfl, 2)), Γη)
  if p.vtk_output
    mkpath(outdir)
    writevtk(Γη, joinpath(outdir, "array_triangulation"),
             cellfields = ["floater" => floater, "checker" => checker])
    writevtk(Λh, joinpath(outdir, "hinge_triangulation"))
    writevtk(Γκ, joinpath(outdir, "surface_triangulation"),
             cellfields = ["μ₁" => CellField(dz.μ₁, Γκ)])
  end

  result = S.simulate(problem, tconfig)

  # Diagnostics: incident η at the array centre, RMS / max array heave and RMS
  # relative hinge rotation |[[∇η]]| (rad)
  xc = Point(lay.xp₀ + Lp / 2, 0.0, 0.0)
  nout = max(1, round(Int, p.Δt_out / p.Δt))
  if p.vtk_output
    io = open(joinpath(outdir, "diagnostics.csv"), "w")
    println(io, "t,eta_inc_centre,eta_rms,eta_max,hinge_rot_rms")
    pvdη = createpvd(joinpath(outdir, "solution_η"))
    pvdκ = createpvd(joinpath(outdir, "solution_κ"))
    pvdϕ = createpvd(joinpath(outdir, "solution_ϕ"))
  end
  for (t, x) in result.solution
    istep = round(Int, t / p.Δt)
    println("solution at t = ", t)
    if p.vtk_output && istep % nout == 0
      ϕₕ, κₕ, ηₕ = x
      η_rms = sqrt(sum(∫(ηₕ * ηₕ)dΓη) / Aarr)
      η_max = maximum(abs, get_free_dof_values(ηₕ))
      θ_rms = sqrt(sum(∫(jump(∇(ηₕ)) ⋅ jump(∇(ηₕ)))dΛh) / Lh)
      println(io, join((t, _ramp(t) * wave.ηin(xc, t), η_rms, η_max, θ_rms), ","))
      flush(io)
      pvdη[t] = createvtk(Γη, joinpath(outdir, "solution_η_$istep"),
                          cellfields = ["η" => ηₕ, "floater" => floater, "checker" => checker]; order = order)
      pvdκ[t] = createvtk(Γκ, joinpath(outdir, "solution_κ_$istep"),
                          cellfields = ["κ" => κₕ, "ηin" => CellField(ηin(t), Γκ)]; order = order)
      pvdϕ[t] = createvtk(Ω, joinpath(outdir, "solution_ϕ_$istep"), cellfields = ["ϕ" => ϕₕ])
    end
  end
  if p.vtk_output
    close(io)
    savepvd(pvdη)
    savepvd(pvdκ)
    savepvd(pvdϕ)
  end

  return nothing
end

function _print_setup(ws, lay, p)
  ℓc = (D / (ρ * g))^(1 / 4)
  println("── Floater array ──────────────────────────────────")
  println("  $nfx × $nfy floaters of $a × $a m  →  $Lp × $Bp m;  ms = $ms kg/m², ",
          "D = $D N·m (hb,eq = $(round(hb, digits=3)) m), ℓc = $(round(ℓc, digits=1)) m, kᵣ = $kᵣ")
  println("── Sea state ──────────────────────────────────────")
  println("  JONSWAP Hs = $Hs m, Tp = $Tp s, γ = $γJ, σθ = $(round(rad2deg(σθ), digits=1))°, H = $H m")
  println("  λp = $(round(ws.λp, digits=1)) m, λmax = $(round(ws.λmax, digits=1)) m, ",
          "λmin = $(round(ws.λmin, digits=1)) m")
  println("── Tank layout ────────────────────────────────────")
  println("  L × W × H = $(lay.L) × $(lay.W) × $H m;  Ld_in = $(lay.Ld_in), Ld_out = $(lay.Ld_out), ",
          "Ld_lat = $(lay.Ld_lat), G_up = $(lay.G_up), G_down = $(lay.G_down), G_lat = $(lay.G_lat)")
  println("── Discretisation ─────────────────────────────────")
  println("  nx × ny × nz = $(lay.nx) × $(lay.ny) × $(p.nz),  Δc = $(p.Δc) m ",
          "($(round(ws.λmin / p.Δc, digits=1)) Q2 cells per λmin), Δf = $(p.Δf) m ",
          "($(round(a / p.Δf, digits=1)) cells per floater),  Δt = $(p.Δt) s ",
          "(Tmin/Δt = $(round(ws.Tmin / p.Δt, digits=1)))")
  ws.λmin / p.Δc < 5 && @warn "Fewer than 5 Q2 cells per shortest wavelength in the coarse zone."
  a / p.Δf < 2 && @warn "Fewer than 2 cells per floater: no C/DG facets inside floaters."
end

# Warm-up (compilation), development run, production run
@timeit to "warm_up" main(WARM)
@timeit to "dev"     main(DEV)
# @timeit to "prod"    main(PROD)

show(to)

end
