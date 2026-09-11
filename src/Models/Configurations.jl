struct PlateParameters{T<:Real}
    fluid_density::T
    gravity::T
    depth::T
    mass_density::T
    tension::T
    bending_rigidity::T
end

struct ResonatorParameters{T<:Real}
    n0::T
    resonator_mass::T
    natural_frequency::T
end

struct FrequencyGradingParameters{T<:Real}
    natural_frequency_start::T
    natural_frequency_minimum::T
    grading_length::T
end

struct DensityGradingParameters{T<:Real}
    n0::T
    grading_length::T
end

"Parameters matching Liu et al. (2025), Table 1."
function liu_2025_parameters(; tension=0.0)
    ρ = 1000.0
    g = 9.8
    depth = 10.0
    β = 0.05
    γ = 0.01
    plate = PlateParameters(ρ, g, depth, ρ * γ, tension, ρ * g * β)
    resonator = ResonatorParameters(1.0, 10.0, 10.0)
    return plate, resonator
end