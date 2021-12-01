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
    "SVK", "SVN", "SWE", "TUR", "TWN", "USA", "RoW"]

industries = 1:S


row_col_names = repeat(ctry_names_2013, inner=S) .* "__" .* string.(repeat(industries, outer=N))
τ = DataFrame([row_col_names ones(N*S, N*S)], ["reporter"; row_col_names])
τ2 = stack(τ, Not(:reporter))
rename!(τ2, :variable => :partner)

# τ_reporter = reduce(vcat, permutedims.(split.(τ2[:,:reporter], "__")))
# τ_partner = reduce(vcat, permutedims.(split.(τ2[:,:variable], "__")))
# τ = DataFrame([τ_reporter τ_partner τ2.variable τ2.value], [:reporter_ctry, :reporter_ind, :partner_ctry, :partner_ind, :variable, :value])



reporter_ctry = ["POL"]
partner_ctry = ["AUT", "BEL"]
reporter_ID = repeat(reporter_ctry .* "__" .* string.(industries), inner=S*length(partner_ctry))
partner_ID = repeat(repeat(partner_ctry, inner=S) .* "__" .* repeat(string.(industries), outer=length(partner_ctry)), outer=S)
τ_new = DataFrame([reporter_ID partner_ID fill(1.2, length(partner_ID))], [:reporter, :partner, :tariffs])


a = leftjoin(τ2, τ_new, on=[:reporter, :partner]) 
a.tariffs .= ifelse.(ismissing.(a.tariffs), a.value, a.tariffs)

b = leftjoin(τ2, a, on=[:reporter, :partner], makeunique=true) # problem now is that it is not sorted!


# takes too long

# row_col_names = repeat(ctry_names_2013, inner=S) .* "__" .* string.(repeat(industries, outer=N))
# τ = DataFrame([row_col_names ones(N*S, N*S)], ["reporter"; row_col_names])
# τ2 = stack(τ, Not(:reporter))
# rename!(τ2, :variable => :partner)


# reporter_ctry = ["POL"]
# partner_ctry = ["AUT", "BEL"]
# reporter_ID = repeat(reporter_ctry .* "__" .* string.(industries), inner=S*length(partner_ctry))
# partner_ID = repeat(repeat(partner_ctry, inner=S) .* "__" .* repeat(string.(industries), outer=length(partner_ctry)), outer=S)
# τ_new = DataFrame([reporter_ID partner_ID fill(1.2, length(partner_ID))], [:reporter, :partner, :tariffs])

# for i in unique(τ_new.reporter)
#     for j in unique(τ_new.partner)
#         tariff = τ_new.tariffs[findfirst((τ_new.reporter .== i) .& (τ_new.partner .== j))]
#         τ2.value[findfirst((τ2.reporter .== i) .& (τ2.partner .== j))] = tariff
#         τ2.value[findfirst((τ2.reporter .== j) .& (τ2.partner .== i))] = tariff
#     end
# end


