function frequency_graded_table(
    plate::PlateParameters,
    resonator::ResonatorParameters,
    grading::FrequencyGradingParameters,
    normalized_wave_numbers,
    normalized_positions;
    profile=ωᵣ_linear,
)
    scale = sqrt(resonator.n0)
    physical_wave_numbers = max.(abs.(normalized_wave_numbers) .* scale, eps(Float64))
    local_frequencies = profile.(normalized_positions .* grading.grading_length,
        Ref(grading))
    branches = [graded_dispersion(plate, resonator, local_frequency, k)
        for local_frequency in local_frequencies, k in physical_wave_numbers]
    return (k_over_sqrt_n0=normalized_wave_numbers,
        x_over_L=normalized_positions, local_frequencies=local_frequencies,
        lower=getproperty.(branches, :lower), upper=getproperty.(branches, :upper))
end

function frequency_graded_positions(
    grading::FrequencyGradingParameters;
    n=7,
    profile=ωᵣ_linear,
)
    endpoint = profile === ωᵣ_linear ?
        1 - grading.natural_frequency_minimum / grading.natural_frequency_start :
        log(grading.natural_frequency_start / grading.natural_frequency_minimum)
    return collect(range(0.0, endpoint; length=n))
end

function frequency_graded_extent(
    grading::FrequencyGradingParameters,
    profile,
)
    if profile === ωᵣ_linear
        return grading.grading_length
    end
    return grading.grading_length *
        log(grading.natural_frequency_start / grading.natural_frequency_minimum)
end