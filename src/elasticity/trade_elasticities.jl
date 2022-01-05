# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script to estimate Caliendo and Parro (2015) trade elasticities
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData, XLSX, LinearAlgebra, Statistics, CSV

dir = "X:/VIVES/1-Personal/Florian/git/IOTables/src/"
include(dir * "model/transform_data.jl") # Script with functions to import and transform raw data
include(dir * "elasticity/elasticities_functions.jl") # Script with functions to compute statistics for estimating sectoral trade elasticities (method of Caliendo and Parro, 2015)
# include(dir * "elasticity/wto_tariffs_2018.jl") # Script to prepare WTO tariff data from 2018

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

dir = "C:/Users/u0148308/Desktop/raw/" # location of raw data

# Data specification

# WIOD rev. 2013
source = "WIOD"
revision = "2013"
year = 2010 # specified year
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

df_tariffs = DataFrame(XLSX.readtable(dir * "WTO/tariff_matrix.xlsx", "Sheet1")...)
τ_Z = Matrix(convert.(Float64, df_tariffs[:,2:end])) # NS×N
τ_Z = 1.0 .+ τ_Z ./ 100 # NS×N

# -------

# CP statistic is computed using absorption (X), i.e. imports not exports
# therefore, we use the transpose of Z, so we have destination as rows and origin as columns as in τ_Z
M = copy(transpose(Z))

lhs_Z, rhs_Z = elasticity_data(M, τ_Z, "intermediate", N, S)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# export regression data with fixed effects included
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

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
col_names = [["industry", "lhs_Z", "rhs_Z"]; ctry_names_2013; "RoW"]

df_reg = DataFrames.DataFrame([FE_ind lhs_Z rhs_Z FE_ctry], col_names)

# CSV.write(dir * "df_reg.csv", df_reg)
# XLSX.writetable(dir * "df_reg.xlsx", df_reg, overwrite=true) # to see results better

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using GLM

df_reg_Z = subset(df_reg, :lhs_Z => ByRow( x -> !(ismissing(x) || isnothing(x) || isnan(x)) ))
df_reg_Z = df_reg_Z[:, [:industry, :lhs_Z, :rhs_Z]]
transform!(df_reg_Z, :industry => ByRow(string), [:lhs_Z, :rhs_Z] .=> ByRow(Float64), renamecols=false) # columns need to have the right type!

XLSX.writetable(dir * "df_reg_Z.xlsx", df_reg_Z, overwrite=true) # to see results better



ols_Z = GLM.lm(@formula(lhs_Z ~ 0 + rhs_Z), df_reg_Z)

# -----

for i in unique(df_reg_Z.industry)
    gdf = subset(df_reg_Z, :industry => ByRow(x-> x == i))
    ols_Z = GLM.lm(@formula(lhs_Z ~ 0 + rhs_Z), gdf)
    println("Result for industry: $i \n $ols_Z \n")
end