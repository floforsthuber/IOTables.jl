# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script to estimate Caliendo and Parro (2015) trade elasticities
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData, XLSX, LinearAlgebra, Statistics, CSV

dir = "X:/VIVES/1-Personal/Florian/git/IOTables/src/"
include(dir * "model/transform_data.jl") # Script with functions to import and transform raw data
include(dir * "elasticity/elasticities_functions.jl") # Script with functions to compute statistics for estimating sectoral trade elasticities (method of Caliendo and Parro, 2015)
# include(dir * "elasticity/wto_tariffs_2018.jl") # Script to prepare WTO tariff data from 2018

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

# Data specification
dir_raw = "C:/Users/u0148308/data/raw/" # location of raw data

# WIOD rev. 2013
source = "WIOD"
revision = "2013"
year = 2011 # specified year
N = 41 # number of countries 
S = 35 # number of industries

# # WIOD rev. 2016
# source = "WIOD"
# revision = "2016"
# year = 2014 # specified year
# N = 44 # number of countries 
# S = 56 # number of industries

# # OECD rev. 2021
# source = "OECD"
# revision = "2021"
# year = 1995 # specified year
# N = 71 # number of countries 
# S = 45 # number of industries

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_data(dir_raw, source, revision, year, N, S)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------


df_tariffs = DataFrame(XLSX.readtable(dir_raw * "WTO/tariff_matrix_" * string(year) * ".xlsx", "Sheet1")...)
τ_Z = Matrix(convert.(Float64, df_tariffs[:,2:end])) # NS×N
τ_Z = 1.0 .+ τ_Z ./ 100 # NS×N

# -------

M = copy(transpose(Z))
Z = copy(M)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# reduce to EU
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

ctry_names_2013 = ["AUS", "AUT", "BEL", "BGR", "BRA", "CAN", "CHN", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "GBR", 
    "GRC", "HUN", "IDN", "IND", "IRL", "ITA", "JPN", "KOR", "LTU", "LUX", "LVA", "MEX", "MLT", "NLD", "POL", "PRT", "ROM", "RUS", 
    "SVK", "SVN", "SWE", "TUR", "TWN", "USA"]

ctry_EU28 = ["AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GBR", "GRC", "HUN", "IRL", "ITA", "LVA", 
    "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROM", "SVK", "SVN", "ESP", "SWE"]

countries = [ctry_names_2013; "ZoW"]

position_EU28 = [2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 25, 26, 27, 29, 30, 31, 32, 33, 35, 36, 37]
position_other = [1, 5, 6, 7, 19, 20, 23, 24, 28, 34, 38, 39, 40, 41]

# -------------

col_others = zeros(N*S, 1)
col_EU = zeros(N*S, S)

for j in 1:N
    if j in position_EU28
        col_EU += Z[:, 1+(j-1)*S:j*S]
    else
        col_others = hcat(col_others, Z[:, 1+(j-1)*S:j*S])
    end
end

Z_temp = hcat(col_others[:,2:end], col_EU)

row_others = zeros(1, size(Z_temp,2))
row_EU = zeros(S, size(Z_temp,2))

for i in 1:N
    if i in position_EU28
        row_EU += Z_temp[1+(i-1)*S:i*S, :]
    else
        row_others = vcat(row_others, Z_temp[1+(i-1)*S:i*S, :])
    end
end

Z_new = vcat(row_others[2:end,:], row_EU)

# -------------

τ_Z_temp = hcat(τ_Z[:, position_other], τ_Z[:, 2])

row_others = zeros(1, size(τ_Z_temp,2))
row_EU = zeros(S, size(τ_Z_temp,2))

for i in 1:N
    if i in position_EU28
        row_EU += τ_Z_temp[1+(i-1)*S:i*S, :]
    else
        row_others = vcat(row_others, τ_Z_temp[1+(i-1)*S:i*S, :])
    end
end

row_EU = row_EU ./ length(position_EU28)

τ_Z_new = vcat(row_others[2:end,:], row_EU)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------


N = 15


lhs_Z, rhs_Z = elasticity_data(Z_new, τ_Z_new, "intermediate", N, S)

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
ctry_names = ["AUS", "BRA", "CAN", "CHN", "IDN", "IND", "JPN", "KOR", "MEX", "RUS", "TUR", "TWN", "USA"]
col_names = [["industry", "lhs_Z", "rhs_Z"]; ctry_names; "RoW"; "EU"]

df_reg = DataFrames.DataFrame([FE_ind lhs_Z rhs_Z FE_ctry], col_names)

CSV.write(dir_raw * "df_reg_EU.csv", df_reg)
XLSX.writetable(dir_raw * "df_reg_EU.xlsx", df_reg, overwrite=true)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using FixedEffectModels, RegressionTables

df_reg_Z = subset(df_reg, :lhs_Z => ByRow( x -> !(ismissing(x) || isnothing(x) || isnan(x)) ))
df_reg_Z = df_reg_Z[:, [:industry, :lhs_Z, :rhs_Z]]
transform!(df_reg_Z, :industry => ByRow(string), [:lhs_Z, :rhs_Z] .=> ByRow(Float64), renamecols=false) # columns need to have the right type!

# XLSX.writetable(dir_raw * "df_reg_Z_EU.xlsx", df_reg_Z, overwrite=true) # to see results better

rr = FixedEffectModels.reg(df_reg_Z, @formula(lhs_Z ~ 0 + rhs_Z), Vcov.robust(), save=true)

for i in unique(df_reg_Z.industry)
    gdf = subset(df_reg_Z, :industry => ByRow(x-> x == i))
    rr = FixedEffectModels.reg(gdf, @formula(lhs_Z ~ 0 + rhs_Z), Vcov.robust(), save=true)
    # reg = RegressionTables.regtable(rr; renderSettings = asciiOutput(), estimformat="%0.4f") # different format
    println("Result for industry: $i \n $rr \n")
end