# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Baseline model
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData,  XLSX, LinearAlgebra, Statistics, CSV

include("transform_data2.jl") # Script with functions to import and transform raw data
include("price_hat.jl") # Script with function to obtain the price index
include("wage_hat.jl") # Script with function to obtain the wages and gross output
include("tariffs_function.jl") # Script with functions to create τ_hat_Z, τ_hat_F from tariff data
include("head_ries_index.jl") # Script with functions to create bilateral Head-Ries index (symmetric bilateral trade costs)
include("elasticities_functions.jl") # Script with functions to compute statistics for estimating sectoral trade elasticities (method of Caliendo and Parro, 2015)
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

dir = "C:/Users/u0148308/Desktop/raw/" # location of raw data

# Data specification

# WIOD rev. 2013
source = "WIOD"
revision = "2013"
year = 1995 # specified year
N = 41 # number of countries 
S = 35 # number of industries

# # WIOD rev. 2016
# source = "WIOD"
# revision = "2016"
# year = 2014 # specified year
# N = 44 # number of countries 
# S = 56 # number of industries


# -------------------------------------------------------------------------------------------------------------------------------------------------------------

Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_data(dir, source, revision, year, N, S)

# -------

# use random tariff matrix for now (tariffs between 0 and 20 percent)
# τ_Z = 1 .+ rand(0.0:0.01:0.2, N*S, N*S) # NS×N
# τ_F = 1 .+ rand(0.0:0.01:0.2, N*S, N*S) # NS×N

df_tariffs = DataFrame(XLSX.readtable(dir * "MFN/tariff_matrix.xlsx", "Sheet1")...)
τ_Z = Matrix(convert.(Float64, df_tariffs[:,2:end]))
τ_Z = 1 .+ τ_Z ./ 100

τ_F = copy(τ_Z)

lhs_Z, rhs_Z, lhs_F, rhs_F = elasticity_data(Z, F, τ_Z, τ_F, N, S)

# -------

combinations = Int64(sum([n*(n+1)/2 for n in 1:N-2]))

# sector FE
industries = lpad.(1:S, 2, '0') # left pad with leading zeros so we can sort later on
FE_ind = repeat(industries, inner=combinations) # S*combinations×1, sectors are already sorted!

# exporter/importer FE
FE_ctry = zeros(combinations, N)
e = 1
for i in 1:N-2 # i -> j
    for j in i+1:N-1 # j -> n
        for n in j+1:N # n -> i
            FE_ctry[e,i] = 1
            FE_ctry[e,n] = 1
            FE_ctry[e,j] = 1
            e += 1
        end
    end
end

FE_ctry = Int64.(repeat(FE_ctry, S)) # S*combinations×N

# combine data in a DataFrame
ctry_names_2013 = ["AUS", "AUT", "BEL", "BGR", "BRA", "CAN", "CHN", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "GBR", 
    "GRC", "HUN", "IDN", "IND", "IRL", "ITA", "JPN", "KOR", "LTU", "LUX", "LVA", "MEX", "MLT", "NLD", "POL", "PRT", "ROM", "RUS", 
    "SVK", "SVN", "SWE", "TUR", "TWN", "USA"]
col_names = [["industry", "lhs_Z", "rhs_Z", "lhs_F", "rhs_F"]; ctry_names_2013; "RoW"]

df_reg = DataFrames.DataFrame([FE_ind lhs_Z rhs_Z lhs_F rhs_F FE_ctry], col_names)

CSV.write(dir * "df_reg.csv", df_reg)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------


using GLM

df_reg_Z = filter(:lhs_Z => x -> !(ismissing(x) || isnothing(x) || isnan(x)), df_reg)
df_reg_Z =  Float64.(df_reg_Z[:, [:lhs_Z, :rhs_Z]]) # columns need to have the right type! 

ols_Z = GLM.lm(@formula(lhs_Z ~ 0 + rhs_Z), df_reg_Z)

# -------

df_reg_F = filter(:lhs_F => x -> !(ismissing(x) || isnothing(x) || isnan(x)), df_reg)
df_reg_F =  Float64.(df_reg_F[:, [:lhs_F, :rhs_F]]) # columns need to have the right type! 

ols_F = GLM.lm(@formula(lhs_F ~ 0 + rhs_F), df_reg_F)
