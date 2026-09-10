function write_dispersion_csv(path, tables)
    open(path, "w") do io
        println(io, "k_over_sqrt_n0,bare_T0,lrh_lower_T0,lrh_upper_T0,bare_Tpositive,lrh_lower_Tpositive,lrh_upper_Tpositive")
        for i in eachindex(tables.zero_tension.k_over_sqrt_n0)
            println(io, join((tables.zero_tension.k_over_sqrt_n0[i],
                tables.zero_tension.bare[i], tables.zero_tension.lower[i],
                tables.zero_tension.upper[i], tables.pretensioned.bare[i],
                tables.pretensioned.lower[i], tables.pretensioned.upper[i]), ","))
        end
    end
    return path
end