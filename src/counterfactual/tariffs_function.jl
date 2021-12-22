# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script with functions to create τ_hat_Z, τ_hat_F from tariff data
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
    country_sector_tariffs(ctry_one::Vector, ctry_two::Vector, tariffs_one::DataFrame, tariffs_two::DataFrame, industries::Vector, N::Integer, S::Integer)

The function matches industry bilateral tariffs which can then be used to create a matrix with the function 'create_τ_hat'. Notice that there are 
    two country and tariff vector, i.e. asymmetric bilateral trade costs are possible. If the model only allows symmetric bilateral trade costs 
    simply use the same DataFrame for both the arguments 'tariffs_one' and 'tariffs_two'.

# Arguments
- `ctry_one::Vector`: a vector of countries (ISO3 codes) on which the tarrifs of vector 'tariffs_one' are levied on by countries in vector 'ctry_two' .
- `ctry_two::Vector`: a vector of countries (ISO3 codes) on which the tarrifs of vector 'tariffs_two' are levied on by countries in vector 'ctry_one' .
- `tariffs_one::DataFrame`: a DataFrame carrying information on industry tariffs levied on 'ctry_two' by 'ctry_one'.
- `tariffs_two::DataFrame`: a DataFrame carrying information on industry tariffs levied on 'ctry_one' by 'ctry_two'.
- `industries::Vector{String}`: a vector of industries to match tariff and country data.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `τ_new::DataFrame`: a DataFrame containing all affected country-industry bilateral trade costs.

# Examples
```julia-repl
julia> τ_new = country_sector_tariffs(ctry_exiting_EU, ctry_remaining_EU, df_tariffs, df_tariffs, industries, N, S)
```
"""

function country_sector_tariffs(ctry_one::Vector, ctry_two::Vector, tariffs_one::DataFrame, tariffs_two::DataFrame, industries::Vector, N::Integer, S::Integer)

    # uses countries in vector "ctry_one" as exporters which pay tariffs "tariffs_one" to importer countries "ctry_two"
    # as the importing ctry levy the tariffs depending on the product (sector) of the exporters
    # we need to match with the exporter with the sector of the exporter!
    exporter_ctry = ctry_one
    importer_ctry = ctry_two

    exporter_ID = repeat(exporter_ctry, inner=S) .* "__" .* repeat(industries, outer=length(exporter_ctry))
    importer_ID = repeat(importer_ctry, inner=S) .* "__" .* repeat(industries, outer=length(importer_ctry))

    exporter = repeat(exporter_ID, inner=length(importer_ID))
    importer = repeat(importer_ID, outer=length(exporter_ID))
    exporter_sector = reduce(vcat, permutedims.(split.(exporter, "__")))[:, 2] # to match with tariff data (on sectoral level)
    importer_ctry_F = reduce(vcat, permutedims.(split.(importer, "__")))[:, 1] # to match with final demand country

    setup_one = DataFrame([exporter importer exporter_sector importer_ctry_F], [:exporter, :importer, :exporter_sector, :importer_ctry])
    df_one = rename(tariffs_one, :industry_WIOD_2016 => :exporter_sector)

    τ_one = leftjoin(setup_one, df_one, on=:exporter_sector)

    # -------

    # uses countries in vector "ctry_two" as exporters which pay tariffs "tariffs_two" to importer countries "ctry_one"
    exporter_ctry = ctry_two
    importer_ctry = ctry_one

    exporter_ID = repeat(exporter_ctry, inner=S) .* "__" .* repeat(industries, outer=length(exporter_ctry))
    importer_ID = repeat(importer_ctry, inner=S) .* "__" .* repeat(industries, outer=length(importer_ctry))

    exporter = repeat(exporter_ID, inner=length(importer_ID))
    importer = repeat(importer_ID, outer=length(exporter_ID))
    exporter_sector = reduce(vcat, permutedims.(split.(exporter, "__")))[:, 2] # to match with tariff data (on sectoral level)
    importer_ctry_F = reduce(vcat, permutedims.(split.(importer, "__")))[:, 1] # to match with final demand

    setup_two = DataFrame([exporter importer exporter_sector importer_ctry_F], [:exporter, :importer, :exporter_sector, :importer_ctry])
    df_two = rename(tariffs_two, :industry_WIOD_2016 => :exporter_sector)

    τ_two = leftjoin(setup_two, df_two, on=:exporter_sector)

    # -------

    τ_new = vcat(τ_one, τ_two)

    return τ_new

end


# -------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
    create_τ_hat(df::DataFrame, scenario::Symbol, countries::Vector, N::Integer, S::Integer)

The function maps the industry bilateral country-industry tariffs supplied by the function 'country_sector_tariffs' in the correct cells 
    for to produce the NS×NS matrix τ_hat_Z, and NS×N matrix τ_hat_F used for computing a counterfactual. Notice bilateral trade costs can 
    vary between intermediate inputs and final demand.

# Arguments
- `df::DataFrame`: a DataFrame providing the counterfactual in form of bilateral country-industry tariffs.
- `scenario_Z::Symbol`: specifies the tariff scenario for intermediate input trade costs.
- `scenario_F::Symbol`: specifies the tariff scenario for final demand trade costs.
- `countries::Vector`: a vector of all countries (ISO3 codes) used in the IO table.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `τ_hat_Z::Matrix{Float64}`: NS×NS, origin country-industry destination country-industry intermediate input trade cost (changes) matrix τ_hat_Z.
- `τ_hat_F::Matrix{Float64}`: NS×N, origin country-industry destination country final demand trade cost (changes) matrix τ_hat_Z.

# Examples
```julia-repl
julia> τ_hat_Z, τ_hat_F = create_τ_hat(df_affected, :avg_tariff_Z_soft, :avg_tariff_F_soft, ctry_names_2013, N, S)
```
"""


function create_τ_hat(df::DataFrame, scenario_Z::Symbol, scenario_F::Symbol, countries::Vector, N::Integer, S::Integer)

    col_names = repeat(countries, inner=S) .* "__" .* repeat(industries, outer=N) # NS×1

    # τ_hat_Z
    df_Z = df[:, [:exporter, :importer, :exporter_sector, :importer_ctry, scenario_Z]]
    rename!(df_Z, scenario_Z => :tariffs)
    
    τ_hat_Z = DataFrame([col_names ones(N*S, N*S)], ["exporter"; col_names])
    τ_hat_Z = stack(τ_hat_Z, Not(:exporter))
    rename!(τ_hat_Z, :variable => :importer)

    τ_hat_Z = leftjoin(τ_hat_Z, df_Z, on=[:exporter, :importer])
    τ_hat_Z.tariffs = ifelse.(ismissing.(τ_hat_Z.tariffs), τ_hat_Z.value, τ_hat_Z.tariffs)
    sort!(τ_hat_Z, [:exporter, :importer]) # sort to have same sorting as IO Table (thats why we use "ZoR" instead of "RoW")
    τ_hat_Z = τ_hat_Z[:, [:exporter, :importer, :tariffs]]
    τ_hat_Z = unstack(τ_hat_Z, :importer, :tariffs, allowduplicates=true)

    τ_hat_Z = Matrix(convert.(Float64, τ_hat_Z[:, 2:end])) # NS×NS

    # -------

    # τ_hat_F
    df_F = df[:, [:exporter, :importer, :exporter_sector, :importer_ctry, scenario_F]]
    rename!(df_F, scenario_F => :tariffs)

    τ_hat_F = DataFrame([col_names ones(N*S, N)], ["exporter"; countries])
    τ_hat_F = stack(τ_hat_F, Not(:exporter))
    rename!(τ_hat_F, :variable => :importer_ctry)

    τ_hat_F = leftjoin(τ_hat_F, df_F, on=[:exporter, :importer_ctry])
    τ_hat_F.tariffs = ifelse.(ismissing.(τ_hat_F.tariffs), τ_hat_F.value, τ_hat_F.tariffs)
    sort!(τ_hat_F, [:exporter, :importer_ctry]) # sort to have same sorting as IO Table (thats why we use "ZoR" instead of "RoW")
    τ_hat_F = τ_hat_F[:, [:exporter, :importer_ctry, :tariffs]]
    τ_hat_F = unstack(τ_hat_F, :importer_ctry, :tariffs, allowduplicates=true) # for some reason it gives out duplicate error?

    τ_hat_F = Matrix(convert.(Float64, τ_hat_F[:, 2:end])) # NS×N

    return τ_hat_Z, τ_hat_F
end
