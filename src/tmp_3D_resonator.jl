module TMP3D_RESONATOR

# 3D time-domain response of a floating Kirchhoff-Love plate carrying an array
# of point mass-spring-damper resonators (locally resonant metaplate) in a
# short-crested JONSWAP sea.
#
# Same tank, damping zones and incident-sea treatment as `tmp_3D.jl`; the only
# additions are the resonators:
#
# - Geometry: one `G.ResonatorDomain` per resonator on the plate surface `:Γη`.
#   `get_integration_domains` turns each into a `DiracDelta{0}` on the 2D plate
#   manifold embedded in 3D (stored under `:δ_i`).
# - Unknowns: one `ConstantFESpace` (VectorValue{1}) per resonator, as in
#   `HydroElasticFEM.Physics.ResonatorArray`.
# - Equations (resonator i at point pᵢ, q = absolute vertical displacement):
#     plate     … + (Kᵢ/ρ)(η - qᵢ) + (Cᵢ/ρ)(ηₜ - qₜᵢ)   at pᵢ
#     resonator (Mᵢ/ρ) qₜₜᵢ + (Cᵢ/ρ)(qₜᵢ - ηₜ) + (Kᵢ/ρ)(qᵢ - η) = 0
#   The resonator equation is scaled by 1/ρ, like the plate equation, so the
#   coupling blocks are symmetric.

import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import WaveSpec as WS
using Gridap
using Gridap.ODEs
using Gridap.FESpaces: ConstantFESpace
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
D  = E * hb^3 / 12 / (1 - ν^2)
d₀ = ρb * hb / ρ
C  = P.build_KL_tensor(E, ν, hb, ρ)

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

# Graded natural frequency along the wave direction, as in
# tmp_graded_resonators.jl (Liu et al. 2025, §5.1): mass fixed, stiffness
# Kᵢ = M ωrᵢ², ωr rising linearly in x from ωr_lo at the upwave edge of the
# plate to ωr_hi at the downwave edge (all resonators in one column, same x,
# share ωr).  The band is centred on the spectral peak ωp with the same
# relative width as the 2D case (8–11 rad/s ≈ ±16 % around 9.5 rad/s).
ωp = 2π / Tp
ωr_lo = 0.85 * ωp    # 0.534 rad/s (T = 11.8 s), upwave column  (ix = 1)
ωr_hi = 1.15 * ωp    # 0.723 rad/s (T =  8.7 s), downwave column (ix = nrx)

# ─────────────────────────────────────────────────────────────────────────────
# Discretisation presets (see tmp_3D.jl for the rationale)
# ─────────────────────────────────────────────────────────────────────────────
# SMOKE: coarse mesh, one vertical cell, few wave components, 2 time steps.
#        Only checks that the 3D resonator problem builds, assembles and steps.
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

"""
    resonator_array(lay)

Locations and per-resonator (M, K, C, ωr) of the `nrx × nry` resonators, one
at the centre of each `Lp/nrx × Bp/nry` plate patch, graded in x.
"""
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
    # One DiracDelta per resonator on the plate surface Γη
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

function main(p; outdir = joinpath(@__DIR__, "..", "data", "tmp_3D_resonator"))
  ws  = wave_scales(p.fmax_fac)
  lay = tank_layout(ws, p.Δ)
  res = resonator_array(lay)
  dz  = damping_coefficients(lay, ws)
  _print_setup(ws, lay, res, p)

  order = 2
  tank   = build_tank(lay, res, p.nz, p.zbase)
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
  δ = [dom[s] for s in res.δsyms]

  wave = _wave(ws, p)
  vin(t)  = (r = _ramp(t); x -> r * wave.vin(x, t))
  ηd(t)   = (r = _ramp(t); x -> r * dz.μ₂(x) * wave.ηin(x, t))
  ∇ₙϕd(t) = (r = _ramp(t); x -> r * dz.μ₁(x) * wave.vzin(x, t))

  # Spaces: [ϕ, κ, η, q₁, …, q_N]
  reffe = ReferenceFE(lagrangian, Float64, order)
  Vϕ = TestFESpace(Ω, reffe; conformity = :H1)
  Vκ = TestFESpace(Γκ, reffe; conformity = :H1)
  Vη = TestFESpace(Γη, reffe; conformity = :H1)
  Vq = [ConstantFESpace(Γη; field_type = VectorValue{1, Float64}) for _ in 1:res.N]
  Y = MultiFieldFESpace([Vϕ, Vκ, Vη, Vq...])
  X = MultiFieldFESpace([TrialFESpace(V) for V in (Vϕ, Vκ, Vη, Vq...)])

  βₕ = 0.5
  γₜ = 0.5
  βₜ = 0.25
  αₕ = γₜ / (βₜ * p.Δt) / g * (1.0 - βₕ) / βₕ
  γ  = 1.0 * order * (order + 1) / dom[:h_η]

  # Resonator terms. `rel(x, i)` is the relative motion η - qᵢ at pᵢ; testing
  # with (v - ξᵢ) gives the plate and resonator equations at once.
  î = VectorValue(1.0)
  Mρ, Kρ, Cρ = res.M ./ ρ, res.K ./ ρ, res.C ./ ρ
  rel(x, i) = x[3] - x[3 + i] ⋅ î
  m_res(x, y) = sum(Mρ[i] * δ[i]((x[3 + i] ⋅ î) * (y[3 + i] ⋅ î)) for i in 1:res.N)
  c_res(x, y) = sum(Cρ[i] * δ[i](rel(y, i) * rel(x, i)) for i in 1:res.N)
  k_res(x, y) = sum(Kρ[i] * δ[i](rel(y, i) * rel(x, i)) for i in 1:res.N)

  ∇ₙ(ϕ) = ∇(ϕ) ⋅ VectorValue(0.0, 0.0, 1.0)
  m(t, x, y) = ∫( d₀ * x[3] * y[3] )dΓη + m_res(x, y)
  c(t, x, y) =
    ∫( βₕ * x[1] * (y[2] + αₕ * y[1]) - x[2] * y[1] )dΓκ +
    ∫( x[1] * y[3] - x[3] * y[1] )dΓη + c_res(x, y)
  function k(t, x, y)
    ϕ, κ, η = x[1], x[2], x[3]
    w, u, v = y[1], y[2], y[3]
    ∫( ∇(w) ⋅ ∇(ϕ) )dΩ +
    ∫( βₕ * (u + αₕ * w) * (g * κ) - dz.μ₂ * κ * w + dz.μ₁ * ∇ₙ(ϕ) * (u + αₕ * w) )dΓκ +
    ∫( v * (g * η) + ∇∇(v) ⊙ (C ⊙ ∇∇(η)) )dΓη +
    ∫( -jump(∇(v)) ⊙ (mean(C ⊙ ∇∇(η)) ⋅ nΛη.⁺) -
       (mean(C ⊙ ∇∇(v)) ⋅ nΛη.⁺) ⊙ jump(∇(η)) +
       D / ρ * γ * jump(∇(v)) ⊙ jump(∇(η)) )dΛη +
    k_res(x, y)
  end
  function b(t, y)
    w, u = y[1], y[2]
    ∫( w * (vin(t) ⋅ nΓin) )dΓin +
    ∫( w * dz.χgen * (vin(t) ⋅ nΓlat) )dΓlat +
    sum(∫( ∇ₙϕd(t) * (u + αₕ * w) - ηd(t) * w )dΓ for dΓ in dΓgen)
  end

  op = TransientLinearFEOperator((k, c, m), b, X, Y; constant_forms = (true, true, true))
  ode_solver = Newmark(LUSolver(), p.Δt, γₜ, βₜ)

  x₀ = zero(X(0.0))
  v₀ = zero(X(0.0))
  a₀ = zero(X(0.0))
  xₜ = solve(ode_solver, op, 0.0, p.T, (x₀, v₀, a₀))

  # Time series of plate and resonator displacement at every resonator point:
  # δᵢ(ηₕ) is the point value η(pᵢ); qᵢ is the single free DOF of its space.
  mkpath(outdir)
  # Resonator metadata for post-processing (ParaView: scripts/lrmm_render_3D.py)
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
    pvdη = createpvd(joinpath(outdir, "solution_η"))
    pvdκ = createpvd(joinpath(outdir, "solution_κ"))
  end
  for (t, x) in xₜ
    istep = round(Int, t / p.Δt)
    κₕ, ηₕ = x[2], x[3]
    ηp = [sum(δ[i](ηₕ)) for i in 1:res.N]
    qp = [get_free_dof_values(x[3 + i])[1] for i in 1:res.N]
    println(io, t, "  ", join(("$(ηp[i]) $(qp[i])" for i in 1:res.N), "  "))
    println("t = $t: max|η(pᵢ)| = $(maximum(abs, ηp)), max|qᵢ| = $(maximum(abs, qp))")
    if p.vtk_output && istep % nout == 0
      # Resonator displacements qᵢ and plate displacement η(pᵢ) at the
      # attachment points, each as one N-component field (constant over the
      # plate); component i follows the row order of resonators.csv.
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
          "(Tr = $(round(2π/ωr_lo, digits=2)) → $(round(2π/ωr_hi, digits=2)) s; ωp = $(round(ωp, digits=3)) rad/s), ",
          "K = $(round(minimum(res.K), digits=1)) → $(round(maximum(res.K), digits=1)) N/m")
end

# Smoke test only (2 steps on a coarse mesh); DEV/PROD are real simulations.
@timeit to "smoke" main(SMOKE)
@timeit to "dev"  main(DEV)
# @timeit to "prod" main(PROD)

show(to)

end
