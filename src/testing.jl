# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script to estimate Caliendo and Parro (2015) trade elasticities
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData, XLSX, LinearAlgebra, Statistics, CSV

dir = "X:/VIVES/1-Personal/Florian/git/IOTables/src/"
include(dir * "model/transform_data.jl") # Script with functions to import and transform raw data
# include(dir * "elasticity/elasticities_functions.jl") # Script with functions to compute statistics for estimating sectoral trade elasticities (method of Caliendo and Parro, 2015)
# include(dir * "elasticity/wto_tariffs_2018.jl") # Script to prepare WTO tariff data from 2018

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

dir = "C:/Users/u0148308/data/raw/" # location of raw data

# Data specification

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

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_data(dir, source, revision, year, N, S)

# -------

ctry_names = ["AUS", "AUT", "BEL", "BGR", "BRA", "CAN", "CHN", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "GBR", 
    "GRC", "HUN", "IDN", "IND", "IRL", "ITA", "JPN", "KOR", "LTU", "LUX", "LVA", "MEX", "MLT", "NLD", "POL", "PRT", "ROM", "RUS", 
    "SVK", "SVN", "SWE", "TUR", "TWN", "USA"]

ctry_EU = ["AUT", "BEL", "BGR", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GBR", "GRC", "HUN", "IRL", "ITA", "LVA", 
    "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROM", "SVK", "SVN", "ESP", "SWE"] # with GBR but without HVR

countries = [ctry_names; "ZoW"]

industries = lpad.(1:S, 2, '0') # left pad with leading zeros so we can sort later on

# -------


function index_subset(N::Integer, S::Integer, 
    ctry_all::Vector{String}, ctry_subset::Vector{String}, ind_all::Vector{String}, ind_subset::Vector{String})
    
    # subsetting rows
    N_subset = length(ctry_subset)
    S_subset = length(ind_subset)

    all = repeat(ctry_all, inner=S) .* "_" .* repeat(ind_all, outer=N) # all row country-industry pairs
    subset = repeat(ctry_subset, inner=S_subset) .* "_" .* repeat(ind_subset, outer=N_subset) # subset of country-industry pairs
    index_subset = findall(in(subset), all)

    return index_subset
end

function ctry_aggregation(Z::Matrix, F::Matrix, N::Integer, S::Integer, 
    ctry_all::Vector{String}, ctry_subset::Vector{String}, ind_all::Vector{String}, ind_subset::Vector{String})

    # obtain indices for subsetting Z and F
    index_row_subset = index_subset(N, S, ctry_all, ctry_subset, ind_all, ind_subset)
    index_col_subset = index_subset(N, S, ctry_all, ctry_subset, ind_all, ind_subset)
    index_col_subset_F = findall(in(ctry_subset), ctry_all) # only have countries in columns

    # subset columns, construct aggregate and rebuild matrix
    N_col_subset = length(ctry_subset)
    S_col_subset = length(ind_subset)

    Z_col_subset = Z[:, index_col_subset]
    Z_col_subset = [sum(Z_col_subset[i,j:S_col_subset:(N_col_subset-1)*S_col_subset+j]) for i in 1:N*S, j in 1:S_col_subset] # sum over industries to create aggregate
    Z_col_other = Z[:, Not(index_col_subset)] # other countries
    Z_new = hcat(Z_col_other, Z_col_subset) # rebuild Z with new aggregate as last columns

    F_col_subset = F[:, index_col_subset_F]
    F_col_subset = [sum(F_col_subset[i,:]) for i in 1:N*S]
    F_col_other = F[:, Not(index_col_subset_F)] # other countries
    F_new = hcat(F_col_other, F_col_subset) # rebuild F with new aggregate as last column

    # subset rows from new matrix, construct aggregate and rebuild matrix
    N_row_subset = length(ctry_subset)
    S_row_subset = length(ind_subset)
    N_new = size(F_new, 2) # new number of countries

    Z_row_subset = Z_new[index_row_subset, :]
    Z_row_subset = [sum(Z_row_subset[i:S_row_subset:(N_row_subset-1)*S_row_subset+i,j]) for i in 1:S_row_subset, j in 1:N_new*S_row_subset]
    Z_row_other = Z_new[Not(index_row_subset), :]
    Z_new = vcat(Z_row_other, Z_row_subset)

    F_row_subset = F_new[index_row_subset, :]
    F_row_subset = [sum(F_row_subset[i:S_row_subset:(N_row_subset-1)*S_row_subset+i,j]) for i in 1:S_row_subset, j in 1:N_new]
    F_row_other = F_new[Not(index_row_subset), :]
    F_new = vcat(F_row_other, F_row_subset)

    return Z_new, F_new
end


a, b = aggregation(Z, F, N, S, countries, ctry_EU, industries, industries)



# function subset(N::Integer, S::Integer, 
#     ctry_all::Vector{String}, ctry_subset::Vector{String}, ind_all::Vector{String}, ind_subset::Vector{String})
    
#     # subsetting rows
#     N_subset = length(ctry_subset)
#     S_subset = length(ind_subset)

#     all = repeat(ctry_all, inner=S) .* "_" .* repeat(ind_all, outer=N) # all row country-industry pairs
#     subset = repeat(ctry_subset, inner=S_subset) .* "_" .* repeat(ind_subset, outer=N_subset) # subset of country-industry pairs
#     index_subset = findall(in(subset), all)

#     return index_subset
# end

# row_index_EU_Z = subset(N, S, countries, ctry_EU, industries, industries)
# col_index_EU_Z = subset(N, S, countries, ctry_EU, industries, industries)
# col_index_EU_F = findall(in(ctry_EU), countries)

# col_n_subset = length(ctry_EU)
# col_s_subset = length(industries)

# row_n_subset = length(ctry_EU)
# row_s_subset = length(industries)


# Z_col_subset = Z[:, col_index_EU_Z]
# Z_col_EU = [sum(Z_col_subset[i,j:col_s_subset:(col_n_subset-1)*col_s_subset+j]) for i in 1:N*S, j in 1:col_s_subset]

# Z_col_NOT = Z[:, Not(col_index_EU_Z)]

# Z_col_new = hcat(Z_col_NOT, Z_col_EU) # now last 35 cols are EU aggregate

# n_NOT = length(countries[findall(!in(ctry_EU), countries)])

# Z_row_subset = Z_col_new[row_index_EU_Z,:]
# Z_row_EU = [sum(Z_row_subset[i:row_s_subset:(row_n_subset-1)*row_s_subset+i,j]) for i in 1:row_s_subset, j in 1:(n_NOT+1)*row_s_subset] # only EU to EU intermediate

# Z_row_NOT = Z_col_new[Not(row_index_EU_Z),:]

# Z_new = vcat(Z_row_NOT, Z_row_EU)

