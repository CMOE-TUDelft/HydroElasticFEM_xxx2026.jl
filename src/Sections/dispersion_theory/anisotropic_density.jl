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

function anisotropic_density_table(
    plate::PlateParameters,
    resonator::ResonatorParameters,
    normalized_wave_numbers,
    normalized_positions;
    Lx=4.0,
    Ly=1.0,
    axis=:x,
)
    axis in (:x, :y) || throw(ArgumentError("axis must be :x or :y"))
    Lx > 0 || throw(DomainError(Lx, "Lx must be positive"))
    Ly > 0 || throw(DomainError(Ly, "Ly must be positive"))
    scale = sqrt(resonator.n0)
    physical_wave_numbers = max.(abs.(normalized_wave_numbers) .* scale,
        eps(Float64))
    length_scale = axis === :x ? Lx : Ly
    densities = [anisotropic_density(
        axis === :x ? position * length_scale : 0.0,
        axis === :y ? position * length_scale : 0.0,
        resonator.n0, Lx, Ly) for position in normalized_positions]
    branches = [density_graded_dispersion(plate, resonator, density, k)
        for density in densities, k in physical_wave_numbers]
    return (k_over_sqrt_n0=normalized_wave_numbers,
        normalized_positions=normalized_positions, densities=densities,
        lower=getproperty.(branches, :lower), upper=getproperty.(branches, :upper),
        axis=axis, Lx=Lx, Ly=Ly)
end

function anisotropic_dispersion_field(
    plate::PlateParameters,
    resonator::ResonatorParameters,
    normalized_wave_numbers;
    Lx=4.0,
    Ly=1.0,
    x_values=collect(range(-4.0, 4.0; length=161)),
    y_values=collect(range(-4.0, 4.0; length=161)),
)
    Lx > 0 || throw(DomainError(Lx, "Lx must be positive"))
    Ly > 0 || throw(DomainError(Ly, "Ly must be positive"))
    scale = sqrt(resonator.n0)
    wave_numbers = max.(abs.(normalized_wave_numbers) .* scale, eps(Float64))
    lower = Array{Float64}(undef, length(y_values), length(x_values),
        length(wave_numbers))
    upper = similar(lower)
    reference = [density_graded_dispersion(plate, resonator, 0.0, k)
        for k in wave_numbers]
    for (row, y) in enumerate(y_values), (column, x) in enumerate(x_values)
        density = anisotropic_density(x * Lx, y * Ly, resonator.n0, Lx, Ly)
        branches = [density_graded_dispersion(plate, resonator, density, k)
            for k in wave_numbers]
        lower[row, column, :] = getproperty.(branches, :lower)
        upper[row, column, :] = getproperty.(branches, :upper)
    end
    lower_reference = reshape(getproperty.(reference, :lower), 1, 1, :)
    upper_reference = reshape(getproperty.(reference, :upper), 1, 1, :)
    return (x=x_values, y=y_values, k_over_sqrt_n0=normalized_wave_numbers,
        lower=lower, upper=upper,
        lower_shift=lower .- lower_reference,
        upper_shift=upper .- upper_reference, Lx=Lx, Ly=Ly)
end