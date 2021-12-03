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

# iteration parameters
vfactor = 0.4
tolerance = 1e-6
max_iteration = 100
iteration = 0
max_error = 1e7

# sectoral trade elasticity
θ = 5 # assumption
θ = fill(θ, N*S) # NS×1, work with country-industry elasticities

# trade costs
τ_hat_Z = ones(N*S, N*S) # NS×NS
τ_hat_F = ones(N*S, N) # NS×N

# initialize wages and price indices
w_hat = ones(N) # N×1

# # adjust trade balance, (if active => adjustments such that there is no trade deficit, i.e. counterfactual in itself)
# TB_ctry .= 0.0


# -------------------------------------------------------------------------------------------------------------------------------------------------------------


# countries present in WIOD
ctry_names_2013 = ["AUS", "AUT", "BEL", "BGR", "BRA", "CAN", "CHN", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "GBR", 
    "GRC", "HUN", "IDN", "IND", "IRL", "ITA", "JPN", "KOR", "LTU", "LUX", "LVA", "MEX", "MLT", "NLD", "POL", "PRT", "ROM", "RUS", 
    "SVK", "SVN", "SWE", "TUR", "TWN", "USA"]
ctry_names_2013 = [ctry_names_2013; "ZoW"] # use "ZoW" instead of "RoW" so we can sort later on

industries = lpad.(1:S, 2, '0') # left pad with leading zeros so we can sort later on

ctry_names_2016 = ["AUS", "AUT", "BEL", "BGR", "BRA", "CAN", "CHE", "CHN", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", 
    "GBR", "GRC", "HRV", "HUN", "IDN", "IND", "IRL", "ITA", "JPN", "KOR", "LTU", "LUX", "LVA", "MEX", "MLT", "NLD", "NOR", "POL", 
    "PRT", "ROU", "RUS", "SVK", "SVN", "SWE", "TUR", "TWN", "USA"]
ctry_names_2016 = [ctry_names_2016; "ZoW"] # use "ZoW" instead of "RoW" so we can sort later on

ctry_EU27 = ["AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", 
    "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE"]


# Tariff data (matched to WIOD rev. 2016 sectoral aggregation)
df_tariffs = DataFrame(XLSX.readtable(dir * "EU_avg_tariffs.xlsx", "Sheet1")...)
transform!(df_tariffs, :industry_WIOD_2016 => ByRow(x->lpad(x, 2, '0')), renamecols=false)
df_tariffs[:, :avg_tariff_F_soft] .= (df_tariffs.avg_tariff .+ 2.77) ./ 100 .+ 1.0
df_tariffs[:, :avg_tariff_F_hard] .= (df_tariffs.avg_tariff .+ 8.31) ./ 100 .+ 1.0
df_tariffs[:, :avg_tariff_Z_soft] .= (df_tariffs.avg_tariff_F_soft .- 1) ./ 2 .+ 1
df_tariffs[:, :avg_tariff_Z_hard] .= (df_tariffs.avg_tariff_F_hard .- 1) ./ 2 .+ 1


# Functions to insert new trade costs into the appropriate cells
function create_country_sector_pairs(party_one::Vector, party_two::Vector, industries::Vector, N::Integer, S::Integer)

    # uses countries in vector "party_one" as exporters
    exporter_ctry = party_one
    importer_ctry = party_two

    exporter_ID = repeat(exporter_ctry, inner=S) .* "__" .* repeat(industries, outer=length(exporter_ctry))
    importer_ID = repeat(importer_ctry, inner=S) .* "__" .* repeat(industries, outer=length(importer_ctry))

    exporter = repeat(exporter_ID, inner=length(importer_ID))
    importer = repeat(importer_ID, outer=length(exporter_ID))
    exporter_sector = reduce(vcat, permutedims.(split.(exporter, "__")))[:, 2] # to match with tariff data (on sectoral level)
    importer_ctry_F = reduce(vcat, permutedims.(split.(importer, "__")))[:, 1] # to match with final demand country

    τ_new_1 = DataFrame([exporter importer exporter_sector importer_ctry_F], [:exporter, :importer, :exporter_sector, :importer_ctry])

    # -------

    # uses countries in vector "party_two" as exporters
    exporter_ctry = party_two
    importer_ctry = party_one

    exporter_ID = repeat(exporter_ctry, inner=S) .* "__" .* repeat(industries, outer=length(exporter_ctry))
    importer_ID = repeat(importer_ctry, inner=S) .* "__" .* repeat(industries, outer=length(importer_ctry))

    exporter = repeat(exporter_ID, inner=length(importer_ID))
    importer = repeat(importer_ID, outer=length(exporter_ID))
    exporter_sector = reduce(vcat, permutedims.(split.(exporter, "__")))[:, 2] # to match with tariff data (on sectoral level)
    importer_ctry_F = reduce(vcat, permutedims.(split.(importer, "__")))[:, 1] # to match with final demand

    τ_new_2 = DataFrame([exporter importer exporter_sector importer_ctry_F], [:exporter, :importer, :exporter_sector, :importer_ctry])

    # -------

    # allows for asymmetric bilateral trade costs (however, model does not!)
    τ_new = vcat(τ_new_1, τ_new_2)

    return τ_new

end


function create_τ_hat(affected::DataFrame, tariffs::DataFrame, scenario::Symbol, countries::Vector, N::Integer, S::Integer)

    # match tariff data with affected
    df = rename(tariffs, :industry_WIOD_2016 => :exporter_sector)
    df = leftjoin(affected, df, on=:exporter_sector)
    df = df[:, [:exporter, :importer, :exporter_sector, :importer_ctry, scenario]]
    rename!(df, scenario => :tariffs)

    col_names = repeat(countries, inner=S) .* "__" .* repeat(industries, outer=N) # NS×1
    
    # τ_hat_Z
    τ_hat_Z = DataFrame([col_names ones(N*S, N*S)], ["exporter"; col_names])
    τ_hat_Z = stack(τ_hat_Z, Not(:exporter))
    rename!(τ_hat_Z, :variable => :importer)

    τ_hat_Z = leftjoin(τ_hat_Z, df, on=[:exporter, :importer])
    τ_hat_Z.tariffs = ifelse.(ismissing.(τ_hat_Z.tariffs), τ_hat_Z.value, τ_hat_Z.tariffs)
    sort!(τ_hat_Z, [:exporter, :importer]) # sort to have same sorting as IO Table (thats why we use "ZoR" instead of "RoW")
    τ_hat_Z = τ_hat_Z[:, [:exporter, :importer, :tariffs]]
    τ_hat_Z = unstack(τ_hat_Z, :importer, :tariffs)

    τ_hat_Z = Matrix(convert.(Float64, τ_hat_Z[:, 2:end])) # NS×NS

    # τ_hat_F
    τ_hat_F = DataFrame([col_names ones(N*S, N)], ["exporter"; countries])
    τ_hat_F = stack(τ_hat_F, Not(:exporter))
    rename!(τ_hat_F, :variable => :importer_ctry)

    τ_hat_F = leftjoin(τ_hat_F, df, on=[:exporter, :importer_ctry])
    τ_hat_F.tariffs = ifelse.(ismissing.(τ_hat_F.tariffs), τ_hat_F.value, τ_hat_F.tariffs)
    τ_hat_F = sort(τ_hat_F, [:exporter, :importer_ctry]) # sort to have same sorting as IO Table (thats why we use "ZoR" instead of "RoW")
    τ_hat_F = τ_hat_F[:, [:exporter, :importer_ctry, :tariffs]]
    τ_hat_F = unstack(τ_hat_F, :importer_ctry, :tariffs, allowduplicates=true) # for some reason it gives out duplicate error?

    τ_hat_F = Matrix(convert.(Float64, τ_hat_F[:, 2:end])) # NS×N

    return τ_hat_Z, τ_hat_F
end

ctry_exiting_EU = ["DEU"]
ctry_remaining_EU = filter(x-> !(x in ctry_exiting_EU), ctry_EU27)

df_affected = create_country_sector_pairs(ctry_exiting_EU, ctry_remaining_EU, industries, N, S)

τ_hat_Z, τ_hat_F = create_τ_hat(df_affected, df_tariffs, :avg_tariff_F_soft, ctry_names_2013, N, S)

XLSX.writetable("clean/tau_Z.xlsx", Tables.table(τ_hat_Z), overwrite=true)
XLSX.writetable("clean/tau_F.xlsx", Tables.table(τ_hat_F), overwrite=true)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------


while max_error > tolerance && iteration <= max_iteration

    # Price indices
    P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat = create_price_index_hat(w_hat, τ_hat_Z, τ_hat_F, VA_coeff, γ, π_Z, π_F, θ)

    # Counterfactual trade shares
    π_hat_Z = ifelse.(isinf.(π_hat_Z), 0.0, π_hat_Z) # remove Inf
    π_hat_F = ifelse.(isinf.(π_hat_F), 0.0, π_hat_F)

    global π_prime_Z = π_Z .* π_hat_Z # NS×NS
    global π_prime_F = π_F .* π_hat_F # NS×N

    # Labor market clearing
    w_hat_prev = copy(w_hat) # store last wage in case new optimization obtains negative wages

    w_hat, Y_prime, Z_prime, F_prime, ETB_ctry = create_wages_hat(w_hat, vfactor, π_prime_Z, π_prime_F, VA_ctry, TB_ctry, γ, α)

    # update iteration parameters
    error = abs.(w_hat .- w_hat_prev)
    max_error = maximum(error) # update error
    iteration += 1 # update iteration count

    if minimum(w_hat) < 0.0
        w_hat[:] = copy(w_hat_prev)
        println("Outer loop: Iteration $iteration completed with error $max_error (wage negative, rerun with previous estimate)")
    else
        println("Outer loop: Iteration $iteration completed with error $max_error")
    end


end


# Country consumer price index and real wage
P0_F_ctry = (ones(S,N) ./ α).^α # S×N
P0_F_ctry = [prod(P0_F_ctry[:,j]) for j in 1:N] # N×1, initial consumer price index

P_F_ctry = (P_hat_F ./ α).^α # S×N, counterfactual consumer price index
P_F_ctry = [prod(P_F_ctry[:,j]) for j in 1:N] # N×1

P_hat_F_ctry = P_F_ctry ./ P0_F_ctry # N×1, changes in country price index

w_hat_real = w_hat ./ P_hat_F_ctry # N×1, real wage changes

# -------------------------------------------------------------------------------------------------------------------------------------

# Gross output
Y_ctry = [sum(Y[i:i+S-1]) for i in 1:S:N*S] # N×1
Y_prime_ctry = [sum(Y_prime[i:i+S-1]) for i in 1:S:N*S] # N×1

Y_hat_ctry = Y_prime_ctry ./ Y_ctry # N×1, nominal

# how to obtain real change since Y is made out of final and intermediate goods?
Y_hat_ctry_real = Y_hat_ctry ./ P_hat_F_ctry