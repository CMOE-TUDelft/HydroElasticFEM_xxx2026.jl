function dispersion_table(
    plate::PlateParameters,
    resonator::ResonatorParameters,
    normalized_wave_numbers;
    include_resonator=true,
)
    scale = sqrt(resonator.n0)
    physical_wave_numbers = max.(abs.(normalized_wave_numbers) .* scale, eps(Float64))
    bare = bare_dispersion.(Ref(plate), physical_wave_numbers)
    if !include_resonator
        return (k_over_sqrt_n0=normalized_wave_numbers, bare=bare,
            lower=nothing, upper=nothing)
    end
    branches = metaplate_dispersion.(Ref(plate), Ref(resonator), physical_wave_numbers)
    return (k_over_sqrt_n0=normalized_wave_numbers, bare=bare,
        lower=getproperty.(branches, :lower), upper=getproperty.(branches, :upper))
end

function bandgap_limits(
    plate::PlateParameters,
    resonator::ResonatorParameters,
    normalized_wave_numbers;
    y_max=15.0,
)
    table = dispersion_table(plate, resonator, normalized_wave_numbers)
    gap_low = maximum(table.lower)
    gap_high = minimum(table.upper)
    margin = 0.10 * (gap_high - gap_low)
    return (max(0.0, gap_low - margin), min(y_max, gap_high + margin))
end