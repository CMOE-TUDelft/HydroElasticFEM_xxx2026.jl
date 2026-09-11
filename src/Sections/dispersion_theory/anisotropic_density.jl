function anisotropic_density(
    x::Real,
    y::Real,
    n0::Real,
    Lx::Real,
    Ly::Real;
    profile=:gaussian,
)
    Lx > 0 || throw(DomainError(Lx, "Lx must be positive"))
    Ly > 0 || throw(DomainError(Ly, "Ly must be positive"))
    ρ = sqrt((x / Lx)^2 + (y / Ly)^2)
    if profile === :gaussian
        return n0 * exp(-ρ^2 / 2)
    end
    throw(ArgumentError("unsupported anisotropic density profile: $profile"))
end

function anisotropic_resonant_wavenumber(
    x::Real,
    y::Real,
    plate::PlateParameters,
    resonator::ResonatorParameters,
    Lx::Real,
    Ly::Real;
    damping_ratio=0.05,
)
    damping_ratio > 0 || throw(DomainError(damping_ratio, "damping ratio must be positive"))
    density = anisotropic_density(x, y, resonator.n0, Lx, Ly)
    mass_loading = density * resonator.resonator_mass
    ωᵣ = resonator.natural_frequency
    effective_mass = plate.mass_density + mass_loading * (1 + im / (2 * damping_ratio))
    coefficient = plate.gravity * plate.fluid_density - ωᵣ^2 * effective_mass
    if plate.bending_rigidity == 0
        K = ωᵣ^2 / (plate.depth * coefficient)
    else
        throw(ArgumentError("anisotropic resonant closed form requires zero bending rigidity"))
    end
    if plate.tension > 0
        K = (-coefficient + sqrt(coefficient^2 +
            4 * plate.tension * ωᵣ^2 / (plate.fluid_density * plate.depth))) /
            (2 * plate.tension / plate.fluid_density)
    end
    wavenumber = sqrt(K)
    return imag(wavenumber) >= 0 ? wavenumber : -wavenumber
end

function anisotropic_density_field(
    plate::PlateParameters,
    resonator::ResonatorParameters;
    Lx=4.0,
    Ly=1.0,
    damping_ratio=0.05,
    x_values=collect(range(-4Lx, 4Lx; length=161)),
    y_values=collect(range(-4Ly, 4Ly; length=161)),
)
    density = [anisotropic_density(x, y, resonator.n0, Lx, Ly)
        for y in y_values, x in x_values]
    attenuation = [imag(anisotropic_resonant_wavenumber(x, y, plate, resonator,
        Lx, Ly; damping_ratio=damping_ratio)) for y in y_values, x in x_values]
    return (x=x_values, y=y_values, density=density, attenuation=attenuation,
        Lx=Lx, Ly=Ly, damping_ratio=damping_ratio)
end

function anisotropic_optical_depths(
    plate::PlateParameters,
    resonator::ResonatorParameters;
    Lx=4.0,
    Ly=1.0,
    damping_ratio=0.05,
    n_points=2001,
)
    x = collect(range(-4Lx, 4Lx; length=n_points))
    y = collect(range(-4Ly, 4Ly; length=n_points))
    attenuation_x = [imag(anisotropic_resonant_wavenumber(
        coordinate, 0.0, plate, resonator, Lx, Ly;
        damping_ratio=damping_ratio)) for coordinate in x]
    attenuation_y = [imag(anisotropic_resonant_wavenumber(
        0.0, coordinate, plate, resonator, Lx, Ly;
        damping_ratio=damping_ratio)) for coordinate in y]
    integrate(values, coordinates) = sum((values[1:end-1] .+ values[2:end]) .*
        diff(coordinates) ./ 2)
    return (A_x=integrate(attenuation_x, x), A_y=integrate(attenuation_y, y),
        x=x, y=y, Lx=Lx, Ly=Ly, attenuation_x=attenuation_x,
        attenuation_y=attenuation_y)
end