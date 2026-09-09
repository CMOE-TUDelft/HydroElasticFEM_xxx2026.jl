module HydroElasticFEM_xxx2026

export PlateParameters, ResonatorParameters, effective_mass,
    bare_dispersion, metaplate_dispersion, liu_2025_parameters

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

"Return the finite-depth fluid added mass per unit area at wavenumber k."
function fluid_added_mass(plate::PlateParameters, k::Real)
    k > 0 || throw(DomainError(k, "wavenumber must be positive"))
    return plate.fluid_density / (k * tanh(k * plate.depth))
end

"Return the plate restoring coefficient K(k)."
function restoring_coefficient(plate::PlateParameters, k::Real)
    k >= 0 || throw(DomainError(k, "wavenumber must be non-negative"))
    return plate.gravity * plate.fluid_density +
        plate.tension * k^2 + plate.bending_rigidity * k^4
end

"Return the resonator dynamic mass per unit area from the appendix model."
function effective_mass(resonator::ResonatorParameters, omega::Real)
    ratio = omega / resonator.natural_frequency
    return resonator.n0 * resonator.resonator_mass / (1 - ratio^2)
end

"Explicit bare-plate frequency omega(k), including finite-depth fluid loading."
function bare_dispersion(plate::PlateParameters, k::Real)
    mass = plate.mass_density + fluid_added_mass(plate, k)
    return sqrt(restoring_coefficient(plate, k) / mass)
end

"Return the two explicit appendix-A metaplate branches at wavenumber k."
function metaplate_dispersion(
    plate::PlateParameters,
    resonator::ResonatorParameters,
    k::Real,
)
    mass = plate.mass_density + fluid_added_mass(plate, k)
    stiffness = restoring_coefficient(plate, k)
    omega_r_squared = resonator.natural_frequency^2
    total_resonator_mass = resonator.n0 * resonator.resonator_mass
    coefficient = stiffness + omega_r_squared * (mass + total_resonator_mass)
    discriminant = coefficient^2 - 4 * omega_r_squared * stiffness * mass
    discriminant >= 0 || throw(DomainError(discriminant, "branches are not real"))
    root = sqrt(max(discriminant, zero(discriminant)))
    lower_squared = (coefficient - root) / (2 * mass)
    upper_squared = (coefficient + root) / (2 * mass)
    return (lower=sqrt(max(lower_squared, zero(lower_squared))),
        upper=sqrt(max(upper_squared, zero(upper_squared))))
end

"Parameters matching the dimensional values reported in Liu et al. (2025), Table 1."
function liu_2025_parameters(; tension=0.0)
    fluid_density = 1000.0
    gravity = 9.8
    depth = 10.0
    beta = 0.05
    gamma = 0.01
    plate = PlateParameters(fluid_density, gravity, depth,
        fluid_density * gamma, tension, fluid_density * gravity * beta)
    resonator = ResonatorParameters(1.0, 10.0, 10.0)
    return plate, resonator
end

end