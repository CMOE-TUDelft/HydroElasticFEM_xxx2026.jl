"""
Baseline parameter sets of the verification test matrix (§0.1), in the
density-normalised notation of the paper: every structural quantity is divided
by the fluid density ρ (subscript ρ), e.g. m_ρ = m/ρ, T_ρ = T/ρ, C_ρ = D/ρ.

`dim = 2` is the vertical section (1D structure, quantities per unit width,
resonators are spanwise rows); `dim = 3` is a 3D fluid with a 2D plate.
"""
module ParameterSets

using Parameters

export HostSpec, ResonatorSpec, ParameterSet
export resonator_from_frequency, natural_frequency, damping_ratio, with_damping_ratio
export resonator_density, P1, P1LF, P2, P3, P4

"""
Bare host: pre-tensioned Kirchhoff–Love plate on a fluid of depth `H`
(`H = Inf` for deep water). `τ`, `τ_b` are the tension and bending damping
times; `ν` is only used by the 3D FE plate.
"""
@with_kw struct HostSpec
    g::Float64 = 9.8
    ρ::Float64 = 1000.0
    H::Float64 = 10.0
    m_ρ::Float64
    C_ρ::Float64 = 0.0
    T_ρ::Float64 = 0.0
    τ::Float64 = 0.0
    τ_b::Float64 = 0.0
    ν::Float64 = 0.0
end

"Single resonator: mass, stiffness and damping per unit fluid density."
@with_kw struct ResonatorSpec
    m_ρ::Float64
    k_ρ::Float64
    c_ρ::Float64 = 0.0
end

"Resonator with natural frequency `ω_r` and damping ratio `ζ`, eq. (resonator_natural_freq_damping)."
function resonator_from_frequency(m_ρ::Real, ω_r::Real; ζ::Real=0.0)
    k_ρ = m_ρ * ω_r^2
    return ResonatorSpec(m_ρ=m_ρ, k_ρ=k_ρ, c_ρ=2ζ * sqrt(k_ρ * m_ρ))
end

natural_frequency(r::ResonatorSpec) = sqrt(r.k_ρ / r.m_ρ)
damping_ratio(r::ResonatorSpec) = r.c_ρ / (2sqrt(r.k_ρ * r.m_ρ))
with_damping_ratio(r::ResonatorSpec, ζ::Real) = resonator_from_frequency(r.m_ρ, natural_frequency(r); ζ=ζ)

"""
A named baseline set. `a` is the resonator spacing (line array, `dim = 2`) or
the square-lattice constant (`dim = 3`); `NaN` for single-resonator sets.
`extras` holds geometry and case data that are not host/resonator constants.
"""
@with_kw struct ParameterSet
    name::String
    dim::Int
    host::HostSpec
    resonator::Union{Nothing,ResonatorSpec}
    a::Float64 = NaN
    extras::Dict{String,Any} = Dict{String,Any}()
end

"Resonators per unit length (`dim = 2`) or per unit area (`dim = 3`, square lattice)."
resonator_density(ps::ParameterSet) = ps.dim == 2 ? 1 / ps.a : 1 / ps.a^2

"P1: Table 1 of the paper (`tab:dispersion_test_parameters`), ω_r = 10 rad/s, a = 1 m."
function P1(; T_ρ::Real=0.0, dim::Int=2, ζ::Real=0.0)
    host = HostSpec(g=9.8, ρ=1000.0, H=10.0, m_ρ=0.01, C_ρ=0.49, T_ρ=T_ρ)
    return ParameterSet(name="P1", dim=dim, host=host,
                        resonator=resonator_from_frequency(0.01, 10.0; ζ=ζ), a=1.0)
end

"P1-LF: low-frequency variant for bathymetry tests, ω_r = 1 rad/s, a = 5 m (n₀m_ρ,r = 0.01)."
function P1LF(; T_ρ::Real=0.0, H::Real=10.0, ζ::Real=0.0)
    host = HostSpec(g=9.8, ρ=1000.0, H=H, m_ρ=0.01, C_ρ=0.49, T_ρ=T_ρ)
    return ParameterSet(name="P1LF", dim=2, host=host,
                        resonator=resonator_from_frequency(0.05, 1.0; ζ=ζ), a=5.0)
end

"""
P2: 50 m floating membrane with one resonator (Results section draft):
LDPE, h = 0.05 m, ρ_m = 910 kg/m³, T = 10 kN/m, ρ = 1025 kg/m³, H = 25 m;
resonator M = 113.75 kg, K = 28.4375 N/m per unit width (ω_r = 0.5 rad/s).
g is not recorded in the draft; the HydroElasticFEM default 9.81 is used.
"""
function P2()
    ρ = 1025.0
    host = HostSpec(g=9.81, ρ=ρ, H=25.0, m_ρ=910.0 * 0.05 / ρ, C_ρ=0.0, T_ρ=1.0e4 / ρ)
    resonator = ResonatorSpec(m_ρ=113.75 / ρ, k_ρ=28.4375 / ρ)
    extras = Dict{String,Any}("L_m" => 50.0, "x_plate" => (50.0, 100.0),
                              "L_domain" => 150.0, "x_1" => 75.0)
    return ParameterSet(name="P2", dim=2, host=host, resonator=resonator, extras=extras)
end

"P3: frozen VT3D-B channel set (`params_VT3D_B.toml`): Table 1 plate, W = 2 m channel, one point resonator."
function P3()
    g = 9.80665
    host = HostSpec(g=g, ρ=1000.0, H=10.0, m_ρ=0.01, C_ρ=0.05g, T_ρ=0.0, ν=0.0)
    extras = Dict{String,Any}("W" => 2.0, "zeta_S4" => 0.02, "omega_S4_over_omega_r" => 0.9,
                              "hb" => 0.01, "rhob" => 1000.0)
    return ParameterSet(name="P3", dim=3, host=host,
                        resonator=resonator_from_frequency(0.01, 10.0), extras=extras)
end

"P4: frozen VT3D-A Yago VLFS set (`params_VT3D_A.toml`), bare plate."
function P4()
    ρ, ρb, hb, E, ν = 1025.0, 256.25, 2.0, 1.19e10, 0.13
    D = E * hb^3 / (12 * (1 - ν^2))
    host = HostSpec(g=9.81, ρ=ρ, H=58.5, m_ρ=ρb * hb / ρ, C_ρ=D / ρ, T_ρ=0.0, ν=ν)
    extras = Dict{String,Any}("L" => 300.0, "B" => 60.0, "hb" => hb, "E" => E, "rhob" => ρb)
    return ParameterSet(name="P4", dim=3, host=host, resonator=nothing, extras=extras)
end

end # module ParameterSets
