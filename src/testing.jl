# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Baseline model
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData,  XLSX, LinearAlgebra, Statistics

include("transform_data.jl") # Script with functions to import and transform raw data
include("price_hat.jl") # Script with function to obtain the price index
include("wage_hat.jl") # Script with function to obtain the wages and gross output

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

years = 2000:2014 # vector with years covered by WIOD Rev. 2016
N = 41 # number of countries 
S = 35 # number of industries
dir = "C:/Users/u0148308/Desktop/raw/" # location of raw data

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_data(dir, "2013", 1995, N, S)


# sectoral trade elasticity
θ = 5 # assumption
θ = fill(θ, N*S) # NS×1, work with country-industry elasticities

# trade costs
τ_hat_Z = ones(N*S, N*S) # NS×NS
τ_hat_F = ones(N*S, N) # NS×N


ctry_names_2013 = ["AUS", "AUT", "BEL", "BGR", "BRA", "CAN", "CHN", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "GBR", 
    "GRC", "HUN", "IDN", "IND", "IRL", "ITA", "JPN", "KOR", "LTU", "LUX", "LVA", "MEX", "MLT", "NLD", "POL", "PRT", "ROM", "RUS", 
    "SVK", "SVN", "SWE", "TUR", "TWN", "USA"]
ctry_names_2013 = [ctry_names_2013; "ZoW"] # use "ZoW" instead of "RoW" so we can sort later on

industries = lpad.(1:S, 2, '0') # left pad with leading zeros so we can sort later on


# tau


function collect_trade_costs(party_one::Vector, party_two::Vector)

    exporter_ctry = party_one
    importer_ctry = party_two

    exporter_ID = repeat(exporter_ctry, inner=S) .* "__" .* repeat(industries, outer=length(exporter_ctry))
    importer_ID = repeat(importer_ctry, inner=S) .* "__" .* repeat(industries, outer=length(importer_ctry))

    exporter = repeat(exporter_ID, inner=length(importer_ID))
    importer = repeat(importer_ID, outer=length(exporter_ID))
    exporter_sector = reduce(vcat, permutedims.(split.(exporter, "__")))[:, 2] # to match with tariff data (on sectoral level)
    importer_ctry_F = reduce(vcat, permutedims.(split.(importer, "__")))[:, 1] # to match with final demand

    τ_new_1 = DataFrame([exporter importer exporter_sector importer_ctry_F fill(1.2, length(exporter_sector)) fill(1.2, length(exporter_sector))], 
        [:exporter, :importer, :exporter_sector, :importer_ctry, :tariffs_Z, :tariffs_F])

    # -------

    exporter_ctry = party_two
    importer_ctry = party_one

    exporter_ID = repeat(exporter_ctry, inner=S) .* "__" .* repeat(industries, outer=length(exporter_ctry))
    importer_ID = repeat(importer_ctry, inner=S) .* "__" .* repeat(industries, outer=length(importer_ctry))

    exporter = repeat(exporter_ID, inner=length(importer_ID))
    importer = repeat(importer_ID, outer=length(exporter_ID))
    exporter_sector = reduce(vcat, permutedims.(split.(exporter, "__")))[:, 2] # to match with tariff data (on sectoral level)
    importer_ctry_F = reduce(vcat, permutedims.(split.(importer, "__")))[:, 1] # to match with final demand

    τ_new_2 = DataFrame([exporter importer exporter_sector importer_ctry_F fill(1.2, length(exporter_sector)) fill(1.2, length(exporter_sector))], 
        [:exporter, :importer, :exporter_sector, :importer_ctry, :tariffs_Z, :tariffs_F])

    # -------

    τ_new = vcat(τ_new_1, τ_new_2)

    return τ_new

end

τ_new = collect_trade_costs(["POL", "GBR"], ["AUT", "BEL"])


# τ_hat_Z

col_names = repeat(ctry_names_2013, inner=S) .* "__" .* repeat(industries, outer=N)
τ = DataFrame([row_col_names ones(N*S, N*S)], ["exporter"; col_names])
τ2 = stack(τ, Not(:exporter))
rename!(τ2, :variable => :importer)

τ3 = leftjoin(τ2, τ_new, on=[:exporter, :importer])
τ3.tariffs_Z = ifelse.(ismissing.(τ3.tariffs_Z), τ3.value, τ3.tariffs_Z)
τ4 = sort(τ3, [:exporter, :importer])
τ4 = τ4[:, [:exporter, :importer, :tariffs_Z]]
τ5 = unstack(τ4, :importer, :tariffs_Z)

τ_hat_Z_new = Matrix(convert.(Float64, τ5[:,2:end])) # NS×NS

# τ_hat_F

col_names = repeat(ctry_names_2013, inner=S) .* "__" .* repeat(industries, outer=N)
τ = DataFrame([row_col_names ones(N*S, N)], ["exporter"; ctry_names_2013])
τ2 = stack(τ, Not(:exporter))
rename!(τ2, :variable => :importer_ctry)

τ3 = leftjoin(τ2, τ_new, on=[:exporter, :importer_ctry])
τ3.tariffs_F = ifelse.(ismissing.(τ3.tariffs_F), τ3.value, τ3.tariffs_F)
τ4 = sort(τ3, [:exporter, :importer_ctry])
τ4 = τ4[:, [:exporter, :importer_ctry, :tariffs_F]]

τ5 = unstack(τ4, :importer_ctry, :tariffs_F, allowduplicates=true)
τ_hat_F_new = Matrix(convert.(Float64, τ5[:,2:end])) # NS×N




