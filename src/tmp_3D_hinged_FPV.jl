module TMP3DHingedFPV

# 3D time-domain hydroelastic response of an array of hinged floaters (offshore
# floating solar, FPV) in a short-crested JONSWAP sea.  Adapted from
# `tmp_3D.jl`; the configuration follows the Oceans of Energy "Nymphaea
# Aurora" system at Hollandse Kust Noord (HKN), Dutch North Sea.
#
# Published data (oceansofenergy.blue/hollandse-kust-noord-1-project, press
# releases June–Oct 2025):
#   0.5 MWp, ≈ 1 ha, 1 400 PV panels on 196 floaters, water depth 25 m,
#   18.5 km off IJmuiden/Egmond aan Zee; survived Storm Amy (9 Bft,
#   Hmax ≈ 6.84 m).  Earlier OoE units: design Hs = 3.5 m, Hmax = 6.5 m.
# Assumed here (not published):
#   196 = 14 × 14 square floaters of a = 7 m (14 · 7 = 98 m ⇒ 0.96 ha),
#   ≈ 7 panels per floater, all neighbours connected by line hinges.
#
# Top view (x = main wave direction, free surface at z = 0):
#
#   y = +W/2 ┌──────┬──────────── lateral zone (+y) ────────────┬────────┐
#            │      │           ┌─┬─┬─┬─┬─┐                     │        │
#            │inlet │  G_up     ├─┼─┼─┼─┼─┤    G_down           │ outlet │
#            │ zone │ ────────▶ ├─┼─┼─┼─┼─┤ ─────────▶          │  zone  │
#            │      │           ├─┼─┼─┼─┼─┤  hinged floaters    │        │
#            │      │           └─┴─┴─┴─┴─┘  (nfx × nfy, a × a) │        │
#   y = -W/2 └──────┴──────────── lateral zone (-y) ────────────┴────────┘
#            0    Ld_in                                        xd_out     L
#
# Structural model
# - The array is one `StructureDomain` carrying a single H¹-continuous η field,
#   so vertical displacement is continuous across every connection line (ideal
#   pinned line hinge).
# - Each floater is an equivalent Kirchhoff-Love plate.  The C/DG slope-
#   continuity terms act only on facets *inside* a floater; the facets on the
#   connection lines are split off into Λh, where the slope may jump freely
#   (zero moment) or through an optional rotational spring kᵣ.
# - Floaters are modelled with the full waterplane (VLFS idealisation):
#   hydrostatic stiffness ρg and areal mass ρ d₀.
#
# Mesh: graded in plan through a piecewise-linear coordinate map, with cell size
# Δf over the array (plus one floater of margin) and Δc elsewhere.  All zone,
# array and hinge lines fall on cell faces.

import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import WaveSpec as WS
using Gridap
using Gridap.ODEs
using TimerOutputs

to = TimerOutput()

# ─────────────────────────────────────────────────────────────────────────────
# Site
# ─────────────────────────────────────────────────────────────────────────────
H  = 25.0      # Water depth at HKN [m]
ρ  = 1025.0    # Sea water density [kg/m³]
g  = 9.81      # Gravity [m/s²]

# ─────────────────────────────────────────────────────────────────────────────
# Floater array (Nymphaea Aurora-type)
# ─────────────────────────────────────────────────────────────────────────────
a   = 7.0      # Floater side length [m]
nfx = 14       # Floaters along x (wave direction)
nfy = 14       # Floaters along y
Lp  = nfx * a  # Array length [m]
Bp  = nfy * a  # Array width  [m]

# Equivalent floater properties.  Areal mass: PV panels (≈ 25 kg each, 7 per
# floater) plus frame, floats and ballast ≈ 100 kg/m².  Bending stiffness of
# the floater frame: D = 1e8 N·m, which makes the characteristic length
# (D/ρg)^¼ ≈ 10 m > a, so each floater is nearly rigid and the global
# flexibility comes from the hinges.
ms  = 100.0                              # Areal mass [kg/m²]
D   = 1.0e8                              # Bending stiffness [N·m]
E   = 70.0e9                             # Aluminium [Pa]
ν   = 0.3
hb  = cbrt(12 * (1 - ν^2) * D / E)       # Equivalent solid thickness [m]
d₀  = ms / ρ                             # Mass per unit area / ρ [m]
C   = P.build_KL_tensor(E, ν, hb, ρ)
kᵣ  = 0.0      # Hinge rotational stiffness per unit hinge length [N·m/rad/m]; 0 = free hinge

# ─────────────────────────────────────────────────────────────────────────────
# Sea state: large but operational wind sea at HKN
# ─────────────────────────────────────────────────────────────────────────────
# Hs = 2.5 m is exceeded roughly 10 % of the year in the southern North Sea and
# stays below the OoE design value (Hs = 3.5 m) and the Storm Amy record.
# Tp = 7.5 s gives a steepness Hs/λp,0 ≈ 1/35, typical of a developed North Sea
# wind sea and still within the range of linear potential-flow theory.
Hs = 2.5       # Significant wave height [m]
Tp = 7.5       # Peak period [s]
γJ = 3.3       # JONSWAP peak enhancement (mean JONSWAP value)
θ̄  = 0.0       # Mean direction [rad]
σθ = deg2rad(25.0)             # Directional spreading of a North Sea wind sea
θmin, θmax = -π / 3, π / 3

# Spectral band kept in the simulation, as multiples of fp
fmin_fac = 0.7

# Damping law (Colomés et al. 2023): μ₁ = μ₀ (1 - cos(π/2 s)), μ₂ = μ₁ kₚ
μ₀ = 2.5

# ─────────────────────────────────────────────────────────────────────────────
# Discretisation presets
# ─────────────────────────────────────────────────────────────────────────────
#   Δc, Δf   coarse / fine (array) horizontal cell size.  a, Lp, Bp must be
#            multiples of Δf; Q2 needs ≥ 5 cells per λmin in the coarse zone.
#   nz,zbase vertical cells and grading base of G.f_z
#   fmax_fac upper spectral cut-off (× fp), sets λmin = λ(fmax)
#   nf, nθ   frequency / direction edges (bins = n - 1)
#   Δt       Newmark step: ≤ Tmin/25 (dev) … Tmin/40 (production)
#   T        simulated time: transient + statistics window
#   Δt_out   VTK / diagnostics output interval
const WARM = (; Δc = 14.0, Δf = 7.0, nz = 1, zbase = 2.0, fmax_fac = 1.2, nf = 20, nθ = 7,
               Δt = 0.3, T = 0.6, Δt_out = 0.3, seed = 1234, vtk_output = true)
const DEV = (; Δc = 7.0, Δf = 3.5, nz = 4, zbase = 2.0, fmax_fac = 1.5, nf = 40, nθ = 11,
               Δt = 0.15, T = 300.0, Δt_out = 0.6, seed = 1234, vtk_output = true)
const PROD = (; Δc = 3.5, Δf = 1.75, nz = 5, zbase = 2.0, fmax_fac = 2.2, nf = 60, nθ = 15,
                Δt = 0.08, T = 1200.0, Δt_out = 0.4, seed = 1234, vtk_output = true)

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

# Round a length up to a whole number of cells
_snap(x, Δ) = Δ * ceil(x / Δ - 1e-9)

_is_multiple(x, Δ) = abs(x / Δ - round(x / Δ)) < 1e-9

"""
    _faces(x₀, segments)

Face coordinates of a 1D graded grid starting at `x₀`, built from
`segments = [(length, Δ), …]` with uniform cells of size `Δ` in each segment.
"""
function _faces(x₀, segments)
  xf = [x₀]
  for (len, Δ) in segments
    _is_multiple(len, Δ) || error("Segment length $len is not a multiple of Δ = $Δ.")
    n = round(Int, len / Δ)
    append!(xf, xf[end] .+ Δ .* (1:n))
  end
  xf
end

# Piecewise-linear map from a uniform grid on [xf[1], xf[end]] with
# length(xf) - 1 cells to the graded faces `xf`
function _face_map(xf)
  n = length(xf) - 1
  x₀, span = xf[1], xf[end] - xf[1]
  ξ -> begin
    s = (ξ - x₀) / span * n
    i = clamp(floor(Int, s), 0, n - 1)
    xf[i + 1] + (s - i) * (xf[i + 2] - xf[i + 1])
  end
end

"""
    tank_layout(ws, p)

Zone lengths derived from the sea state (Tp = 7.5 s, H = 25 m gives
λp ≈ 85 m, λmax ≈ 145 m):

- `Ld_in  ≈ λmax`   inlet generation/absorption zone
- `Ld_out ≈ 2 λp`   outlet absorber
- `Ld_lat ≈ λp`     lateral generation/absorption zones
- `G_up = G_down ≈ λp`  free water between zones and array
- `G_lat ≈ λp/2`    free water between array side and lateral zone

The fine-mesh block covers the array plus one floater `a` of margin.
"""
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

"""
    build_tank(lay, nz, zbase)

3D `TankDomain` with free surface at z = 0, plan-graded and vertically graded
mesh, the floater array as a single `StructureDomain` and four
non-overlapping damping zones.
"""
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
    damping_zones = [
      G.DampingZone(L = Ld_in, x₀ = [0.0, -W / 2, 0.0], domain_symbol = :Γd_in),
      G.DampingZone(L = xd_out - Ld_in, W = Ld_lat, x₀ = [Ld_in, -W / 2, 0.0], domain_symbol = :Γd_ym),
      G.DampingZone(L = xd_out - Ld_in, W = Ld_lat, x₀ = [Ld_in, W / 2 - Ld_lat, 0.0], domain_symbol = :Γd_yp),
      G.DampingZone(L = lay.Ld_out, x₀ = [xd_out, -W / 2, 0.0], domain_symbol = :Γd_out),
    ],
  )
end

"""
    split_hinges(Λη, lay)

Split the array skeleton `Λη` into the facets inside floaters (`Λfl`, C/DG
plate terms) and the facets on the interior connection lines (`Λh`, hinges).
"""
function split_hinges(Λη, lay)
  tol = 1e-6 * a
  xh = [lay.xp₀ + i * a for i in 1:nfx-1]
  yh = [lay.yp₀ + j * a for j in 1:nfy-1]
  on_line(c, lines) = any(l -> abs(c - l) < tol, lines)
  is_hinge(xs) = all(x -> on_line(x[1], xh), xs) || all(x -> on_line(x[2], yh), xs)
  mask = map(is_hinge, collect(get_cell_coordinates(Λη)))
  ih, ifl = findall(mask), findall(!, mask)
  nh_expected = round(Int, ((nfx - 1) * Bp + (nfy - 1) * Lp) / _array_Δ(lay))
  length(ih) == nh_expected ||
    error("Found $(length(ih)) hinge facets, expected $nh_expected.")
  Triangulation(Λη, ifl), Triangulation(Λη, ih)
end

# Cell size over the array (fine block)
_array_Δ(lay) = (i = searchsortedfirst(lay.xf, lay.xp₀ - 1e-9); lay.xf[i + 1] - lay.xf[i])

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

function _wave(ws, p)
  Dθ = WS.AngularSpreading.DiscreteAngularSpreading(:normal, θ̄, σθ, θmin, θmax, p.nθ)
  Cspec = WS.ContinuousSpectrums.JONSWAP(Hs, Tp, γJ)
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
  r = WS.realize(sea_state)
  ηin(x, t) = WS.evaluate_η(r, x[1], x[2], t)
  vzin(x, t) = WS.evaluate_w(r, x[1], x[2], x[3], t)
  vin(x, t) = VectorValue(WS.evaluate_u(r, x[1], x[2], x[3], t),
                          WS.evaluate_v(r, x[1], x[2], x[3], t),
                          vzin(x, t))
  (; ηin, vin, vzin)
end

# Smooth start of the forcing over 3 Tp to avoid an impulsive transient
_ramp(t; Tr = 3 * Tp) = t < Tr ? 0.5 * (1 - cos(π * t / Tr)) : 1.0

# ─────────────────────────────────────────────────────────────────────────────
# Solver
# ─────────────────────────────────────────────────────────────────────────────

function main(p; outdir = joinpath(@__DIR__, "..", "data", "tmp_3D_hinged_FPV"))
  ws  = wave_scales(p.fmax_fac)
  lay = tank_layout(ws, p)
  dz  = damping_coefficients(lay, ws)
  _print_setup(ws, lay, p)

  # Geometry
  #   Ω  fluid volume      Γκ free surface incl. damping zones (κ field)
  #   Γη floater array     Λfl skeleton inside floaters, Λh hinge lines
  order = 2
  tank   = build_tank(lay, p.nz, p.zbase)
  model  = G.build_model(tank)
  trians = G.build_triangulations(tank, model)
  dom    = G.get_integration_domains(trians; degree = 2 * order)
  Ω, Γκ, Γη = trians[:Ω], trians[:Γκ], trians[:Γη]
  dΩ, dΓκ, dΓη, dΓin = dom[:dΩ], dom[:dΓκ], dom[:dΓη], dom[:dΓin]
  Λfl, Λh = split_hinges(trians[:Λη], lay)
  dΛfl, nΛfl = Measure(Λfl, 2 * order), get_normal_vector(Λfl)
  dΛh = Measure(Λh, 2 * order)
  nΓin = get_normal_vector(trians[:Γin])
  Γlat = trians[:Γlateral]
  dΓlat, nΓlat = Measure(Γlat, 2 * order), get_normal_vector(Γlat)
  dΓgen = [Measure(trians[s], 2 * order) for s in (:Γd_in, :Γd_ym, :Γd_yp)]

  Aarr = sum(∫(1.0)dΓη)
  Lh   = sum(∫(1.0)dΛh)
  println("  array area = $(round(Aarr, digits=1)) m², hinge length = $(round(Lh, digits=1)) m")

  # Per-cell floater index (i, j) from the cell centroids; `checker` alternates
  # between neighbouring floaters so they can be told apart in the renderings
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

  wave = _wave(ws, p)
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
  γ  = 1.0 * order * (order + 1) / dom[:h_η]   # C/DG penalty, h = Δf

  # Weak form
  ∇ₙ(ϕ) = ∇(ϕ) ⋅ VectorValue(0.0, 0.0, 1.0)
  m(t, (ϕₜₜ, κₜₜ, ηₜₜ), (w, u, v)) = ∫( d₀ * ηₜₜ * v )dΓη
  c(t, (ϕₜ, κₜ, ηₜ), (w, u, v)) =
    ∫( βₕ * ϕₜ * (u + αₕ * w) - κₜ * w )dΓκ +
    ∫( ϕₜ * v - ηₜ * w )dΓη
  k_plate(v, η) =
    ∫( v * (g * η) + ∇∇(v) ⊙ (C ⊙ ∇∇(η)) )dΓη +
    # C/DG slope continuity inside each floater only (hinge facets excluded)
    ∫( -jump(∇(v)) ⊙ (mean(C ⊙ ∇∇(η)) ⋅ nΛfl.⁺) -
       (mean(C ⊙ ∇∇(v)) ⋅ nΛfl.⁺) ⊙ jump(∇(η)) +
       D / ρ * γ * jump(∇(v)) ⊙ jump(∇(η)) )dΛfl
  # Rotational spring across the hinge lines (moment = kᵣ × relative rotation)
  k_hinge(v, η) = kᵣ > 0 ? ∫( kᵣ / ρ * jump(∇(v)) ⋅ jump(∇(η)) )dΛh : ∫( 0.0 * v * η )dΓη
  k(t, (ϕ, κ, η), (w, u, v)) =
    ∫( ∇(w) ⋅ ∇(ϕ) )dΩ +
    # free surface + relaxation (μ₁ = μ₂ = 0 outside the damping zones)
    ∫( βₕ * (u + αₕ * w) * (g * κ) - dz.μ₂ * κ * w + dz.μ₁ * ∇ₙ(ϕ) * (u + αₕ * w) )dΓκ +
    k_plate(v, η) + k_hinge(v, η)
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

  # Diagnostics: incident η at the array centre, RMS / max array heave and RMS
  # relative hinge rotation |[[∇η]]| (rad)
  xc = Point(lay.xp₀ + Lp / 2, 0.0, 0.0)
  nout = max(1, round(Int, p.Δt_out / p.Δt))
  if p.vtk_output
    mkpath(outdir)
    io = open(joinpath(outdir, "diagnostics.csv"), "w")
    println(io, "t,eta_inc_centre,eta_rms,eta_max,hinge_rot_rms")
    pvdη = createpvd(joinpath(outdir, "solution_η"))
    pvdκ = createpvd(joinpath(outdir, "solution_κ"))
    pvdϕ = createpvd(joinpath(outdir, "solution_ϕ"))
  end
  for (t, x) in xₜ
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
                          cellfields = ["κ" => κₕ, "ηin" => CellField(x -> _ramp(t) * wave.ηin(x, t), Γκ)]; order = order)
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
  println("  $nfx × $nfy floaters of $a × $a m  →  $Lp × $Bp m;  ms = $ms kg/m² (d₀ = $(round(d₀, digits=3)) m), ",
          "D = $D N·m (hb,eq = $(round(hb, digits=3)) m), ℓc = $(round(ℓc, digits=1)) m, kᵣ = $kᵣ")
  println("── Sea state ──────────────────────────────────────")
  println("  JONSWAP Hs = $Hs m, Tp = $Tp s, γ = $γJ, σθ = $(round(rad2deg(σθ), digits=1))°, H = $H m")
  println("  λp = $(round(ws.λp, digits=1)) m, λmax = $(round(ws.λmax, digits=1)) m (f = $(round(ws.fmin, digits=4)) Hz), ",
          "λmin = $(round(ws.λmin, digits=1)) m (f = $(round(ws.fmax, digits=4)) Hz)")
  println("── Tank layout ────────────────────────────────────")
  println("  L × W × H = $(lay.L) × $(lay.W) × $H m;  Ld_in = $(lay.Ld_in), Ld_out = $(lay.Ld_out), ",
          "Ld_lat = $(lay.Ld_lat), G_up = $(lay.G_up), G_down = $(lay.G_down), G_lat = $(lay.G_lat)")
  println("  array x ∈ [$(lay.xp₀), $(lay.xp₀ + Lp)], y ∈ [$(lay.yp₀), $(lay.yp₀ + Bp)]")
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
