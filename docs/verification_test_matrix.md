# Verification test matrix — monolithic FE formulation (J2025-LRMM)

Date: 2026-10-12. Status: **proposal** (not pre-registered). Each test below
must be frozen through the normal pre-registration route before a binding run.

Scope: verify `eq:discreteform_freq` (and `eq:weakform_time`) against
(i) exact references of Sec. 3 (`sec:dispersion_relation`) and Sec. 4
(`sec:discrete_lattice`), (ii) geometry-independent identities, and
(iii) local/WKB asymptotics with controlled error, extending to problems
without analytical solution (3D, complex plan shapes, bathymetry).

Equation labels refer to `dispersion_theory.tex`, `sec_discrete_lattice.tex`,
`lrmm_numerical.tex`, `problem_setting.tex` in the Overleaf project.

---

## 0. Shared definitions

### 0.1 Baseline parameter sets

| Set | Use | g | H | m_ρ | C_ρ | T_ρ | Resonator (per resonator) | Spacing / density | ω_r |
|---|---|---|---|---|---|---|---|---|---|
| **P1** Table-1 plate | L0, L1, L2.5, L3 | 9.8 | 10 m | 0.01 m | 0.49 m⁵/s² | 0 / 9.8 m³/s² | m_ρr = 0.01 m² (1D), 0.01 m³ (3D); k_ρr = m_ρr ω_r² | a = 1 m (n₀ = 1 m⁻¹ or m⁻²) | 10 rad/s |
| **P1-LF** low-frequency variant | L4.1–4.3, L4.6 | 9.8 | 10 m (varies) | 0.01 | 0.49 | 0 | m_ρr = 0.05 m², k_ρr = 0.05 | a = 5 m | 1 rad/s |
| **P2** membrane case | L2.3 | 9.8 | 25 m | 0.0444 m | 0 | 9.7561 | M = 113.75 kg/m, K = 28.4375 N/m/m, c = 0 | single, x₁ = 75 m on plate [50,100] | 0.5 rad/s |
| **P3** 3D-B channel | L3.2 | 9.80665 | 10 m | 0.01 | 0.05 g | 0 | M_ρ = 0.01 m³, ζ = 0.02 | single, W = 2 m channel | 10 rad/s |
| **P4** Yago (3D-A) | L3.1b | 9.81 | 58.5 m | ρ_b h_b/ρ (ρ_b=256.25, h_b=2) | E=1.19e10, ν=0.13 | 0 | — | — | — |

P3 and P4 are the frozen sets in `deliverables/06_Numerical_Verification/data/VT3D_{B,A}/params_*.toml`; reuse them unchanged.

### 0.2 Derived reference values (P1, bare host, lossless; computed from `eq:lat_symbol`)

| ω (rad/s) | T_ρ | k₀ (1/m) | λ_p (m) | λ_open (m) | 𝒟₀′(k₀) | γ_r 1D = 1/𝒟₀′ | γ_r 3D = k₀/2𝒟₀′ | ζ_crit 1D | ζ_crit 3D | Q_rad 1D |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 0 | 1.758 | 3.58 | 2.46 | 18.7 | 0.0534 | 0.0469 | 0.053 | 0.047 | 75 |
| 9.646 | 0 | 2.691 | 2.34 | 0.68 | 51.1 | 0.0196 | 0.0264 | 0.010 | 0.014 | 55 |
| 10 | 0 | 2.742 | 2.29 | 0.62 | 53.7 | 0.0186 | 0.0255 | 0.0093 | 0.0128 | 54 |
| 5 | 9.8 | 1.113 | 5.64 | 2.46 | 44.7 | 0.0224 | 0.0125 | 0.022 | 0.0125 | 179 |
| 10 | 9.8 | 1.926 | 3.26 | 0.62 | 78.7 | 0.0127 | 0.0122 | 0.0064 | 0.0061 | 79 |

Near-field scales at ω = 10: flexural-root decay length 1/Im k ≈ 0.37 m (T=0), 0.23 m (T=9.8); fluid evanescent decay H/π = 3.18 m.
ζ_crit uses |κ_r| ≈ k_ρr (ωc ≪ k_ρr); the script must use the exact condition γ_d = γ_r.

P1-LF wavenumbers (bare host, T=0): ω=1: λ_p = 51.6 m (H=10), 32.3 m (H=3), 59.7 m (H=20); ω=2: 15.3 m (H=10), 13.6 m (H=3).

Published band edges (P1, `fig:lat_bands`): [ω_X, ω_Γ] = [9.646, 9.986] (T=0), [9.916, 9.990] (T=9.8); Π_L = 2.86e-3 / 1.94e-3.

### 0.3 Common discretisation and acceptance rules

| Item | Default |
|---|---|
| FE orders | p_φ = 2; κ, η Q2 with C/DG plate (γ = 6); also p = 3 for rate checks |
| Mesh ladder | h_plate = λ_p/{8,16,32}; open water h ≤ λ_open/6.4; vertical nz = {8,16,32} graded towards surface (top cell k·h_z ≤ 0.2) |
| Resonator near field | local refinement h_min = λ_p/{12,24,48} within R = λ_p/4 of each x_i (3D) |
| Truncation | generating/absorbing zones ≥ 4λ, gap ≥ max(λ/2, 3H/π) between structure/array and zone |
| Identity tests (I) | relative residual ≤ 1e-10 (round-off) |
| Discretisation tests (D) | observed rate ≥ theoretical − 0.3 on the two finest levels, and Richardson-extrapolated error ≤ 1e-3 relative |
| Asymptotic tests (A) | FE − local-theory discrepancy E(ε) fitted as Cε^p; p within ±0.3 of the predicted order, C reported |
| Compute | one solver-heavy job at a time under `Loop_Framework/tools/run_limited.sh`, fixed threads, provenance record |

Observable conventions: R, T normalised by the flux `eq:lat_flux`; Bloch q from the phase of the cell-to-cell ratio η(x+a)/η(x) or from a Bloch unit-cell eigenproblem; complex q at real ω when damping is present.

---

## Level 0 — Discrete-operator consistency

### L0.1 Independence of free-surface stabilisation (α, β)
Theory link: `eq:Bilinear_freq`, coercivity argument of `colomes2023monolithic`. Type: I/D.

| Parameter | Values |
|---|---|
| Geometry | P1, 2D, plate on x ∈ [0, 20] m, open water 2 × 20 m each side + zones |
| Resonators | 5, a = 1, centred, c = 0 |
| ω | {2, 5, 9.8, 10.2} |
| β | {0.1, 0.5, 0.9} |
| α | {α_ref/2, α_ref, 2α_ref} (α_ref from colomes2023monolithic) |
| Mesh | middle level only |

Output: ‖η_h(α,β) − η_h(α_ref,0.5)‖/‖η_h‖, |q_i| variation. Pass: ≤ 1e-8 (solver tol).

### L0.2 Limiting constituents
Theory link: verification list item 1 of `subsec:lat_verification`. Type: I.

| Case | Resonator values | Compared against |
|---|---|---|
| Bare host | k_ρr = 1e-10 k_ref and separately m_ρr = 1e-10 m_ref | FE without resonator DOFs |
| Rigid attachment | k_ρr = 1e10 k_ref | FE with point mass m_ρr added to M₁₁ |
| Dry plate | fluid removed; free edges | `eq:eigDry` modes, analytic Euler–Bernoulli free–free + point mass |
| Geometry | P1, plate 20 m, one resonator at x = 7.3 m (off-node and on-node variants) | |

Pass: relative difference ≤ 1e-6 (penalty-limited) and decreasing linearly with the regularising parameter.

### L0.3 Discrete energy balance in time domain (`prop:energy`)
Theory link: `eq:energy`, `eq:energies`. Type: I (lossless) / D (lossy).

| Parameter | Values |
|---|---|
| Domain | closed tank, length 40 m, H = 10 m, walls (u_in = u_out = 0) |
| Plate | P1 on x ∈ [10, 30], T_ρ ∈ {0, 9.8} |
| Resonators | 10, a = 2 m, ω_r = 10 |
| Initial condition | Gaussian η, width 2 m, amplitude 0.01 m, q_i = η(x_i) |
| Damping cases | (a) τ = τ_b = c = 0; (b) ζ = 0.05; (c) τ = 0.01 s, τ_b = 0.01 s |
| Time integration | implicit midpoint (energy-conserving) and Newmark-β(¼,½); Δt = T_r/{20,40,80}; 100 T_r |
| β | {0.25, 0.75} |

Output: E_t(t), ∫P_t dt, drift. Pass (a): |E_t(t) − E_t(0)|/E_t(0) ≤ 1e-10 for midpoint; (b,c): |E_t(t) − E_t(0) + ∫P_t| ≤ O(Δt²).

### L0.4 Power accounting per resonator (frequency domain)
Theory link: `eq:powDissTime`, `eq:lat_alpha` (Im α ≥ 0). Type: I.

| Parameter | Values |
|---|---|
| Geometry | as L0.1 |
| ζ | {0.01, 0.05, 0.2} |
| ω | 50 points in [8, 12] |

Output: ½ω² c_ρ,i |η(x_i) − q_i|² vs. Im(α_i)|η(x_i)|²/ω-form; open-domain balance P_in = P_R + P_T + P_ds + P_dr. Pass: per-resonator identity 1e-12; global balance treated as D (radiation-boundary error), report rate.

---

## Level 1 — Vertical section, configuration A (fully covered surface)

Common setup: plate covers the whole top boundary; incident hydroelastic wave generated in a generating zone; absorbing zone at the far end tuned to k₀ of the plate (not the open-water k). Gap ≥ 3H/π = 9.6 m between array ends and zones.

### L1.1 Bare-host dispersion
Theory link: `eq:lat_symbol`, `eq:Zf_closed_form`, deep/shallow limits `eq:deep_water_poly`, `eq:shallow_water_poly`. Type: D.

| Parameter | Values |
|---|---|
| Plate | P1, T_ρ ∈ {0, 9.8}; dry-check C_ρ ∈ {0.049, 0.49, 4.9} |
| H | {1, 10, 100} m (shallow / reference / deep) |
| ω | 30 log-spaced in [0.3, 15] |
| Domain | channel length 20λ_p, measurement window 10λ_p |
| Mesh | ladder §0.3, p ∈ {2,3} |

Output: k_h(ω) from phase fit; |k_h − k₀|/k₀. Pass: rate 2p (dispersion error), extrapolated ≤ 1e-4.

### L1.2 Flux normalisation
Theory link: `eq:lat_flux`, `eq:lat_gamma_r`. Type: D.

| Parameter | Values |
|---|---|
| Forcing | unit point force f at x = 0, symmetric domain ±10λ_p + zones |
| ω | {2, 5, 10} |
| T_ρ | {0, 9.8} |

Output: work ½ρωγ_r|F|² vs. 2 × far-field flux from `eq:lat_flux`; Im G_h(0) vs. 1/𝒟₀′(k₀) (table §0.2). Pass: D.

### L1.3 Single line resonator
Theory link: `eq:lat_detuning`, `eq:lat_renorm_single`, `eq:lat_absorbed_1d`, `eq:lat_Qrad`. Type: D.

| Parameter | Values |
|---|---|
| Host | P1, T_ρ ∈ {0, 9.8} |
| Resonator | m_ρr = 0.01, ω_r = 10 |
| ζ | {0 (lossless), 0.5ζ_crit, ζ_crit, 2ζ_crit} with ζ_crit from γ_d = γ_r (≈ 0.0093 / 0.0064) |
| ω | 401 points in ω_r(1 ± 5/Q_rad) plus refinement near peak |
| Extra | attachment on-node vs mid-element |

Outputs and references:

| Observable | Reference |
|---|---|
| Peak frequency ω̃_r | root of `eq:lat_renorm_single` with Re G(0) from residue series `eq:lat_green_1d` |
| R, T at peak, c = 0 | R = −1, T = 0 |
| max 𝒜 at ζ_crit | 1/2 |
| Half-power width | 2/Q_rad, Q_rad ≈ 54 / 79 |

### L1.4 Finite periodic array: Foldy, energy, Fabry–Pérot
Theory link: `eq:lat_foldy`, `eq:lat_RT`, `eq:lat_energy_1d`, `eq:lat_fabry_perot`. Type: I (Foldy) + D (FE).

| Parameter | Values |
|---|---|
| Host | P1, T_ρ ∈ {0, 9.8} |
| N | {1, 2, 5, 10, 20, 40} |
| a | 1 m |
| ζ | {0, 0.01, 0.05} |
| ω | 801 points in [8, 12] (fine in [9.6, 10.0]) |
| Reference | semi-analytic Foldy with residue series, M_roots ∈ {1, 5, 20, 80} |

Outputs: |R|, |T| (FE vs Foldy); energy-identity residual; −ln|T_N|/(Na) vs Im q a from `eq:lat_1d_exact`; Fabry–Pérot peak positions N q a = pπ in pass bands.

### L1.5 Discrete stop band vs homogenised anticrossing (headline)
Theory link: `eq:lat_bloch`, `eq:lat_exact_em`, `eq:lat_gamma`, `subsec:lat_long_wave`, `fig:lat_bands`. Type: D.

| Parameter | Values |
|---|---|
| Host / resonators | P1, T_ρ ∈ {0, 9.8}, a = 1, c = 0 |
| Route A (needs Bloch BC) | unit cell [0, a] × [−H, 0], quasi-periodic φ, η (and plate rotation); q a ∈ 101 points in [0, π]; solve for ω (nonlinear eig) or scan ω for real/complex q |
| Route B (no new code) | N = 40 array from L1.4; edges from onset of exponential decay of |T_N| |
| Mesh | cell nx = {16, 32, 64}, nz = {8,16,32} |

| Quantity | Reference T=0 | Reference T=9.8 |
|---|---|---|
| ω_X | 9.646 | 9.916 |
| ω_Γ | 9.986 (≈ ω_r/√(1+Π_L)) | 9.990 |
| Gap Δω/ω_r | 0.034 | 0.0074 |
| Must **not** match | anticrossing 0.163 | 0.138 |

Pass: band edges to 1e-4 relative (D); also show homogenised FE (L1.7) transmits at every ω.

### L1.6 Spacing sweep: single-mode vs multimode law
Theory link: `eq:lat_1d_exact`, `eq:lat_1d_single`. Type: D/A.

| Parameter | Values |
|---|---|
| a | {0.25, 0.5, 1, 2, 4, 8} m (spans 0.37 m flexural and 3.18 m fluid decay lengths) |
| Resonator | m_ρr = 0.01a (fixed mass ratio n₀m_ρr = 0.01), ω_r = 10 |
| T_ρ | 0 |
| Route | as L1.5 (Route A preferred) |
| Reference | `eq:lat_1d_exact` with M ∈ {1, 3, 10, 40} roots |

Output: band edges vs a; error of single-mode law vs a, expected ∝ e^{−a Im k_m}.

### L1.7 Homogenised FE and the lattice-compliance number Π_L
Theory link: `eq:dispersion_relation`, `eq:omega_k_x_closed_form`, `eq:lat_PiL`, `eq:lat_gamma`. Type: I (homogenised FE) / A.

Code variant: resonator field q(x) ∈ V_Γm with density n(x) replacing point sums (cheap to add).

| Parameter | Values |
|---|---|
| a (point model) | {0.5, 1, 2, 3} m, with m_ρr = 0.01a, k_ρr = a → Π_L ≈ a⁴/352.8 = {1.8e-4, 2.8e-3, 4.5e-2, 0.23} |
| T_ρ | {0, 9.8} |
| Homogenised FE | n₀ m_ρr = 0.01, same ω_r |
| ω | 401 points in [8, 12] |

Outputs: (i) homogenised FE k(ω) vs `eq:omega_k_x_closed_form` (I, discretisation only); (ii) ω_Γ,point − ω_r/√(1+Π_L) (A, expected O(Π_L²)); (iii) difference point-vs-homogenised dispersion vs ε_s = d/λ and Π_L.

---

## Level 2 — Vertical section, configuration B (finite structure in open water)

### L2.1 Transfer compliance G_{B,h} and reciprocity
Theory link: `sec:finite_domain`. Type: D (reciprocity residual is a truncation diagnostic).

| Parameter | Values |
|---|---|
| Plate | P1, length L = 20 m, free edges |
| Attachment points | N = 20 (a = 1) plus 5 random points |
| ω | {2, 5, 9.8, 10, 10.2} |
| Truncation | zone length {1, 2, 4, 8}λ_open, gap {0.25, 0.5, 1}λ_open |

Output: ‖G − Gᵀ‖/‖G‖ vs zone length and mesh. Pass: monotone decrease; ≤ 1e-6 at finest.

### L2.2 Woodbury equivalence
Theory link: `eq:lat_foldy` with G = G_{B,h}. Type: I.

| Parameter | Values |
|---|---|
| Geometry | as L2.1 (single mesh) |
| Resonator sets | (a) identical; (b) random ω_r ∈ [8,12], ζ ∈ [0, 0.1], m_ρr ∈ [0.005, 0.02], seed fixed; (c) linear grading |
| ω | 21 points in [8, 12] |

Pass: ‖η_a,mono − η_a,Foldy‖/‖η_a‖ ≤ 1e-10.

### L2.3 Renormalised wet resonance (existing membrane case)
Theory link: `eq:lat_renorm_single` with G_B(x₁,x₁); modal problem `eq:eigWet`. Type: D.

| Parameter | Values |
|---|---|
| Set | P2 (L_m = 50 m, domain 150 m, H = 25 m) |
| Mesh | dx = {0.8, 0.4, 0.2} m, ny = {8, 15, 30} |
| Observables | complex wet eigenvalue near 0.4997 − 1.6e-5 i; forced-response peak |
| Reference | ω̃ from `eq:lat_renorm_single` with Re G_B(x₁,x₁;ω) and width from Im G_B |
| Extra | x₁ ∈ {50 (edge), 62.5, 75} m |

Pass: eigenvalue, forced peak and Green-function prediction agree to D tolerance; imaginary part vs γ_r = Im G_B.

### L2.4 Absorption bound with asymmetric radiation
Theory link: `eq:lat_absorbed_B`, `eq:lat_absorbed_1d`. Type: D.

| Parameter | Values |
|---|---|
| Plate | P1, L ∈ {5, 20} m |
| x₁ | {centre, L/4, L/8, 0.25 m from edge} |
| ω | ω_r = 10 (and ω̃_r from L2.3-type renormalisation) |
| ζ | sweep 0.001–0.1 (40 log pts) to locate γ_d = γ_r |
| Incidence | from −∞ and from +∞ |

Outputs: a^± from unit-force radiation; max 𝒜 vs bound |a⁻|²/(|a⁻|²+|a⁺|²). Pass: attainment within D tolerance; centre case = 1/2.

### L2.5 Finite-size cut-off of the fluid inertia (μ_B)
Theory link: `eq:lat_mu`, μ_B ~ n₀m_ρr π tanh(πH/L)/L (stated "to be verified"). Type: A.

| Parameter | Values |
|---|---|
| Plate | P1, T_ρ = 0 |
| L / λ_p(ω_r) | {1, 2, 5, 10, 20} → L ≈ {2.3, 4.6, 11.5, 23, 46} m |
| Resonators | a = 1, N = L |
| H | {3, 10, 30} m |
| ω | 401 in [8, 12] |

Output: hybridisation splitting of the coupled wet modes (`eq:eigWet` with resonators) vs μ_B scaling.

### L2.6 Graded arrays (frequency and density grading)
Theory link: `eq:x_star_linear`, `eq:x_star_exp`, `eq:footprint_lin_exp`, `eq:density_grade`, `eq:homogenization_scale_separation`. Type: A. Use Woodbury (L2.2) for parameter sweeps.

| Parameter | Linear | Exponential | Density |
|---|---|---|---|
| Law | ω_r = 10(1 − x/L_g) | ω_r = 10 e^{−x/L_g} | n = n₀e^{−x/L_g}, ω_r = 10; implement as m_ρr,i = m_ρr e^{−x_i/L_g} at fixed a |
| L_g | 500 m | 500 m (array length 1151 m to reach ω_r = 1) | {λ, 5λ, 20λ, 100λ} at ω = 9.8 |
| a | 1 m and 0.25 m (m_ρr scaled, ε_s check) | same | 1 m |
| ζ | {0.01, 0.02, 0.05} | same | {0, 0.02} |
| ω | {2, 3, 5, 7} | {2, 3, 5, 7} | 41 in [9.5, 10.5] |
| Predicted x* | {400, 350, 250, 150} m | {805, 602, 347, 178} m | — |
| Predicted Δx (ζ = 0.02) | {4, 6, 10, 14} m | 20 m (ω-independent) | — |
| Abrupt-end comparison | — | — | uniform array truncated sharply, same N_eff |

Outputs: x of maximum P_dr,i(ω) and its FWHM; |R| for tapered vs abrupt arrays. Note: at a = 1 the high-ω end has ε_s ≈ 0.3–0.44, where the local law is not expected to hold; the a = 0.25 runs separate discreteness error from grading error.

---

## Level 3 — 3D fluid with a 2D plate, flat bottom

### L3.1 Bare plate baseline
(a) Theory link `eq:lat_symbol` with oblique wave vector. (b) existing VT3D-A (P4, Yago) — reuse.

| Parameter | Values |
|---|---|
| Domain (a) | periodic in y, width W = λ_p/sinθ-compatible; θ ∈ {0°, 30°, 60°} |
| Plate | P1, T_ρ = 0 |
| ω | {2, 5, 10} |
| Mesh | ladder §0.3, nz = {8, 16} |

Pass: |k_h − k₀| D tolerance; isotropy (θ-independence).

### L3.2 Point resonator regularisation (existing VT3D-B)
Theory link: Remark in `subsec:governing_eqs_time`; footprint 𝓕(kr_c). Type: D. Use set P3 unchanged (h_min/λ_p = 12, 24, 48; R = λ_p/4).

| Added variant | Values |
|---|---|
| C_ρ | {0.05g, 0.005g, 5e-4g} — show η(x_i) bounded, growth ~ C^{-1/2} |
| Footprint radius r_c | {0, λ_p/50, λ_p/20} |
| Attachment location | node / edge / element interior |

### L3.3 Single resonator: cross-sections and optical theorem
Theory link: `eq:lat_green_2d`, `eq:lat_bound`, `eq:lat_farfield`, `eq:lat_optical`. Type: D.

| Parameter | Values |
|---|---|
| Domain | plate-covered square 16λ_p × 16λ_p (≈ 37 m at ω = 10), H = 10, absorbing frame 4λ_p |
| Resonator | P1 3D, m_ρr = 0.01 m³, centred |
| ζ | {0, ζ_crit ≈ 0.0128, 3ζ_crit} |
| ω | 41 points in ω_r(1 ± 3/Q_rad) |
| Far-field ring | r = {4, 5, 6}λ_p, 72 angles |
| Symmetry | quarter domain for radiation problem only |

| Observable | Reference |
|---|---|
| Im G_h(0) | γ_r = k₀/(2𝒟₀′) = 0.0255 |
| σ_abs max | 1/k₀ = 0.365 m at ζ_crit |
| σ_sc max (c=0) | 4/k₀ = 1.46 m |
| Optical theorem | residual ≤ D tolerance |

### L3.4 Doubly periodic lattice band structure
Theory link: `eq:lat_bloch`, `eq:lat_kummer`, `eq:lat_2d_rect`, `eq:lat_PiL` (square). Type: D. Requires Bloch BCs.

| Parameter | Values |
|---|---|
| Lattice | square a = 1 m; rectangular a = 1, b ∈ {0.5, 2} m |
| Resonator | m_ρr = 0.01 m³, ω_r = 10, c = 0 → Π_L ≈ 3.87e-3/0.49 = 7.9e-3 |
| T_ρ | {0, 9.8} |
| q sampling | path Γ–X–M–Γ (60 pts) **and** full zone 41 × 41 |
| Mesh | cell 16² / 32² in-plane, nz = {8, 16} |
| Reference | `eq:lat_kummer` with β = k₀, |G| ≤ 200 · 2π/a; cross-check `eq:lat_2d_rect` |

Output: ω_Γ, ω_X, ω_M, complete stop band (full-zone), comparison with ω_r/√(1+Π_L).

### L3.5 Grating stack
Theory link: `eq:lat_energy_grating`. Type: I (semi-analytic) / D (FE).

| Parameter | Values |
|---|---|
| Rows N | {1, 3, 10} |
| Row period b | {0.5, 1} m; row spacing a = 1 m |
| θ | {0°, 20°, 45°} |
| ω | 201 in [8, 12] (open orders change at Wood anomalies) |

### L3.6 Finite cluster
Theory link: `eq:lat_foldy` with `eq:lat_green_2d`. Type: D.

| Parameter | Values |
|---|---|
| Cluster | disc of radius 3a, a = 1 (≈ 29 resonators) and random (seed fixed) same count |
| Domain | as L3.3 |
| ω | {9.6, 9.8, 10.0} |
| Reference | Foldy with residue Green function (M_roots = 20) |

Output: f(θ), σ_sc, σ_abs vs Foldy; optical theorem.

---

## Level 4 — Beyond analytical solutions

Verification of every L4 test = (i) identities (energy balance, reciprocity of G_{B,h}, Woodbury, optical theorem/Haskind where applicable), (ii) convergence in h, p, truncation, (iii) agreement with local/WKB theory as ε → 0 with fitted order.

### L4.1 Slowly varying bathymetry, plate-covered
Theory link: `eq:Zf_closed_form` with H(x); `eq:lat_flux`; Sec. 3 statement that uniform n₀ gives x-dependent dispersion. Type: A.

| Parameter | Values |
|---|---|
| Profile | H(x) = H₁ + (H₂ − H₁)[1 + tanh((x − x_c)/L_H)]/2, H₁ = 10, H₂ = 3 m |
| L_H / λ_p | {0.25, 0.5, 1, 2, 4, 8} |
| Plate | P1-LF host (T = 0) |
| Resonators | none (a), then uniform n₀ = 0.2 m⁻¹, ω_r = 1 (b) |
| ω | {0.5, 1, 2} (k₀H spans 0.28–4.1) |
| Dimension | 2D first; 3D oblique incidence θ = 30° second (refraction) |

Outputs: local k(x) from phase gradient vs root of `eq:lat_symbol` at H(x); WKB invariant |A(x)|² 𝒟₀′(k₀; H(x)) = const; reflection |R| → 0 as L_H/λ grows (expected exponential); in 3D, ray angle vs Snell law with local k.

### L4.2 Depth step under resonators
Type: A + identities.

| Parameter | Values |
|---|---|
| Step | H = 10 → 3 m at x = 0 (sharp and L_H = λ_p/4) |
| Resonators | N = 20, a = 5 m, on shallow side, P1-LF (ω_r = 1) |
| ω | 401 in [0.8, 1.2] |
| Reference for band edges | `eq:lat_bloch` evaluated at H = 3 and H = 10 |
| Optional | mode-matching solution for bare step (if available) |

Outputs: stop-band edges vs local-depth Bloch prediction; energy balance; reciprocity.

### L4.3 Graded resonators over graded depth (combined turning point)
Theory link: `eq:turning_point_general` with k(x) from H(x). Type: A.

| Parameter | Values |
|---|---|
| Grading | ω_r = 2 e^{−x/L_g}, L_g = 500 m, a = 5 m |
| Depth | linear ramp 20 → 3 m over [0, 700] m; and constant 10 m (control) |
| ζ | 0.02 |
| ω | {0.6, 0.8, 1.0, 1.4} |

Outputs: x* and footprint vs `eq:x_star_exp`, `eq:footprint_lin_exp`; shift between ramp and constant-depth cases vs local-theory prediction.

### L4.4 Anisotropic Gaussian island (3D-C)
Theory link: `eq:elliptic_density`, `eq:iso_contour_ellipse`, `eq:aspect_ratio_law`. Type: A. Gate: L3.2 and L3.3 pass.

| Parameter | Values |
|---|---|
| Plate / host | P1 (T = 0), H = 10, plate-covered domain with absorbing frame |
| Aspect L_x/L_y | {0.5, 1, 2} at fixed L_xL_y = (6λ_p)² |
| Size scaling | min(L_x, L_y)/λ_p ∈ {3, 6} (tail error expected to drop; report) |
| Sampling | Protocol D, fixed seed, n₀ = 1 m⁻² peak, total resonator mass fixed |
| ω | ω_r = 10 (resonant) and 9.5 |
| ζ | {0.02, 0.05} |
| Incidence | +x and +y |

Outputs: 𝒜_x/𝒜_y from plate field and from resonator-amplitude decay; ellipticity of |η| level sets; convergence of the ratio towards L_x/L_y with island size.

### L4.5 Complex plan shapes
Type: identities + D. No analytical reference.

| Case | Geometry (P1 plate, H = 10, ω ∈ {5, 9.8, 10, 10.2}) | Specific check |
|---|---|---|
| L-shape | three squares of side 4λ_p; re-entrant corner | corner-condition implementation; local refinement convergence near the corner |
| Perforated | square 8λ_p with circular hole radius λ_p | free-edge BC on inner boundary |
| Two plates | two 4λ_p squares, gap {λ_open/4, λ_open} | gap resonance; reciprocity across plates |
| Resonators | 1 resonator scanned over 25 positions; then 25-resonator array | placement map of max 𝒜 / σ_abs ≤ 1/k₀ bound from G_B(x,x) |

Outputs: reciprocity residual of G_{B,h}, Haskind relation (excitation from radiation far field), energy balance; convergence rates.

### L4.6 Annular plate around a monopile (application)
Type: identities + D.

| Parameter | Values |
|---|---|
| Monopile | rigid cylinder radius 3.5 m, bottom-mounted |
| Plate | annulus 3.5 → 50 m, P1 properties, free outer edge, clamped/free inner edge (two variants) |
| H | 20 m (k₀H ≈ 2.1 at ω = 1) |
| Resonators | rings at r = {10, 20, 30, 40} m, angular spacing 5 m; ω_r ∈ {0.8, 1.0, 1.2} |
| ζ | {ζ_crit, 0.05} |
| ω | 41 in [0.5, 1.5] |
| Incidence | plane wave; also directional spread ±30° |

Outputs: load on pile with/without resonators, σ_abs per resonator vs 1/k₀, energy balance, Woodbury-based optimisation reproducibility.

### L4.7 Irregular sea in time domain
Theory link: `eq:weakform_time`, `prop:energy` with radiation. Type: D.

| Parameter | Values |
|---|---|
| Case | L2.6 exponential array (2D) and L4.6 (3D, if budget) |
| Spectrum | JONSWAP, T_p = 2π/ω (ω ∈ {2, 5}), γ = 3.3, H_s = 0.1 m (linear) |
| Duration | 200 T_p, Δt = T_p/80 |
| Absorbing | zones as §0.3 |

Outputs: RAO reconstructed from time series vs frequency-domain FE (coherence ≥ 0.99 in band); cumulative energy budget E_in − E_out − ∫P_t.

---

## Missing capabilities (blockers)

| Capability | Needed by | Note |
|---|---|---|
| Quasi-periodic (Bloch) BCs on φ, κ, η and plate rotation | L1.5-A, L1.6, L3.1a, L3.4, L3.5 | not in `eq:weakform_freq`; Route B fallback for 1D only |
| Absorbing zone tuned to plate-covered k₀ | all config-A tests | current RBC uses one k, exact only for propagating mode |
| Homogenised resonator field | L1.7 | small code variant |
| Rigid bottom-mounted cylinder | L4.6 | wall BC on cylinder surface |
| G_B column assembly + Foldy solver | L2.1–2.6, L4.5–4.6 | one factorisation per ω |
| Semi-analytic references (residue Green function, Kummer lattice sums, Foldy) | L1–L3 | extend `bloch_bands_fig.py` / `src/plate_dispersion.jl` |

## Execution order and existing scripts

1. L0.1–L0.4 (2D, minutes) — extend `src/verification.jl`.
2. L1.1–L1.4 — `case_A_single_resonator.jl`, `case_C_periodic_array.jl`.
3. L2.1–L2.3 — `case_D_finite_structure.jl`, `case_B_two_resonators.jl`.
4. L1.5–L1.7 (needs Bloch BC / homogenised variant).
5. L2.4–L2.6, L4.1–L4.3 (2D) — `case_E_graded_structure.jl`.
6. L3.1–L3.3 — `VT3D_A_yago.jl`, `VT3D_B_point.jl`.
7. L3.4–L3.6, then L4.4 (3D-C), L4.5, L4.6, L4.7 as budget allows.

Each level gates the next: no L(n+1) binding run before L(n) passes its identity tests.
