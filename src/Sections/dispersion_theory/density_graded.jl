function density_graded_positions(; n=7)
    return collect(range(0.0, 1.0; length=n))
end

function density_graded_table(
    plate::PlateParameters,
    resonator::ResonatorParameters,
    grading::DensityGradingParameters,
    normalized_wave_numbers,
    normalized_positions,
)
    scale = sqrt(resonator.n0)
    physical_wave_numbers = max.(abs.(normalized_wave_numbers) .* scale, eps(Float64))
    densities = n_density_graded.(normalized_positions .* grading.grading_length,
        Ref(grading))
    branches = [density_graded_dispersion(plate, resonator, local_density, k)
        for local_density in densities, k in physical_wave_numbers]
    return (k_over_sqrt_n0=normalized_wave_numbers,
        x_over_L=normalized_positions, densities=densities,
        mass_loading=Mᵣ.(Ref(resonator), densities),
        local_frequency=fill(resonator.natural_frequency, length(densities)),
        lower=getproperty.(branches, :lower), upper=getproperty.(branches, :upper))
end