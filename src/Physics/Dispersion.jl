"Return the finite-depth fluid added mass per unit area at wavenumber k."
function fluid_added_mass(plate::PlateParameters, k::Real)
    k > 0 || throw(DomainError(k, "wavenumber must be positive"))
    return plate.fluid_density / (k * tanh(k * plate.depth))
end

"Return the plate restoring coefficient Kᵣ(k)."
function Kᵣ(plate::PlateParameters, k::Real)
    k >= 0 || throw(DomainError(k, "wavenumber must be non-negative"))
    return plate.gravity * plate.fluid_density +
        plate.tension * k^2 + plate.bending_rigidity * k^4
end

const restoring_coefficient = Kᵣ

"Return the resonator dynamic mass per unit area from the appendix model."
function effective_mass(resonator::ResonatorParameters, ω::Real)
    ratio = ω / resonator.natural_frequency
    return resonator.n0 * resonator.resonator_mass / (1 - ratio^2)
end

"Explicit bare-plate frequency ω(k), including finite-depth fluid loading."
function bare_dispersion(plate::PlateParameters, k::Real)
    mass = plate.mass_density + fluid_added_mass(plate, k)
    return sqrt(Kᵣ(plate, k) / mass)
end

"Return the two explicit Appendix-A metaplate branches at wavenumber k."
function metaplate_dispersion(
    plate::PlateParameters,
    resonator::ResonatorParameters,
    k::Real,
)
    mass = plate.mass_density + fluid_added_mass(plate, k)
    stiffness = Kᵣ(plate, k)
    ωᵣ² = resonator.natural_frequency^2
    mᵣ = resonator.n0 * resonator.resonator_mass
    coefficient = stiffness + ωᵣ² * (mass + mᵣ)
    discriminant = coefficient^2 - 4 * ωᵣ² * stiffness * mass
    discriminant >= 0 || throw(DomainError(discriminant, "branches are not real"))
    root = sqrt(max(discriminant, zero(discriminant)))
    lower² = (coefficient - root) / (2 * mass)
    upper² = (coefficient + root) / (2 * mass)
    return (lower=sqrt(max(lower², zero(lower²))),
        upper=sqrt(max(upper², zero(upper²))))
end