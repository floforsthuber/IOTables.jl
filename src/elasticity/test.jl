# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData, XLSX, LinearAlgebra, Statistics, CSV

dir = "X:/VIVES/1-Personal/Florian/git/IOTables/src/"
include(dir * "model/transform_data.jl") # Script with functions to import and transform raw data
include(dir * "model/price_hat.jl") # Script with function to obtain the price index
include(dir * "model/wage_hat.jl") # Script with function to obtain the wages and gross output
include(dir * "counterfactual/tariffs_function.jl") # Script with functions to create τ_hat_Z, τ_hat_F from tariff data
include(dir * "counterfactual/head_ries_index.jl") # Script with functions to create bilateral Head-Ries index (symmetric bilateral trade costs)
include(dir * "elasticity/elasticities_functions.jl") # Script with functions to compute statistics for estimating sectoral trade elasticities (method of Caliendo and Parro, 2015)

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
# prepare tariff data from WTO
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

path = dir * "WTO/EU_2018_MFN_BPT.csv"

df2 = CSV.read(path, DataFrame, delim=',', quoted=true, quotechar='"')

# problem with reading data (one column gets split into multiple)
cols = ["Indicator", "Reporting Economy", "Partner Economy ISO3A Code", "Product/Sector Code", "Year", "Value", "Column25", "Column26", "Column27", "Column28"]
df = df2[:, cols]

transform!(df, :Column28 => ByRow(x -> ifelse(ismissing(x), missing, string(x))), renamecols=false) # tariffs imported as Int
# remove all strings which are split over columns and fill with NaN instead
transform!(df, [:Value, :Column25, :Column26, :Column27, :Column28] .=> ByRow(x -> ifelse(ismissing(x) || !contains(x, " "), x, NaN)), renamecols=false)

# function to extract the tariff value or missing of from the split columns
function extract_tariff_value(a, b, c, d, e)
    v = [a, b, c, d, e]
    index = findall(typeof.(v) .== String)
    value = ifelse(any(typeof.(v) .== String), parse.(Float64, v[index]), [missing])
    return value[1]
end

df.tariff = extract_tariff_value.(df.Value, df.Column25, df.Column26, df.Column27, df.Column28)

# format table
df = df[:, Not(["Value", "Column25", "Column26", "Column27", "Column28"])]
df.Indicator .= ifelse.(df.Indicator .== "HS MFN - Simple average ad valorem duty", "MFN", "BPT")
rename!(df, ["indicator", "reporter_ctry", "partner_ctry", "product_code", "year", "tariff"])
transform!(df, ["indicator", "reporter_ctry", "partner_ctry"] .=> ByRow(string), renamecols=false)
df.partner_ctry .= ifelse.(df.indicator .== "MFN", "WTO_MFN", df.partner_ctry)

# use MFN tariff for tariffs which are missing
df_MFN = subset(df, :indicator => x -> x .== "MFN")
df_BPT = subset(df, :indicator => x -> x .== "BPT")

df_comb = leftjoin(df_BPT, df_MFN[:,[:product_code, :tariff]], on=:product_code, makeunique=true) # possibly reversing join would yield more observations
# but at least we would not lose the tariffs for product aggregates
rename!(df_comb, ["tariff", "tariff_1"] .=> ["BPT", "MFN"])
df_comb.tariff = ifelse.(ismissing.(df_comb.BPT), df_comb.MFN, df_comb.BPT)

df_formated = df_comb[:, Not(["indicator", "BPT", "MFN"])]

# only 3.5% of tariffs are different to 0?!
# possibly better if we restrict ctry list
count(df_formated.tariff .== 0.0)/nrow(df_formated)





# -------------------------------------------------------------------------------------------------------------------------------------------------------------

ctry = ["AUS", "CAN", "CHN", "IDN", "JPN", "KOR", "TWN", "USA"]

function transform_wto_data(path::String)
    
    df = CSV.read(path, DataFrame, delim=',', quoted=true, quotechar='"')
    cols = ["Indicator", "Reporting Economy ISO3A Code", "Partner Economy ISO3A Code", "Product/Sector Code", "Year", "Value"]
    df = df[:, cols]

    df.Indicator .= ifelse.(df.Indicator .== "HS MFN - Simple average ad valorem duty", "MFN", "BPT")
    rename!(df, ["indicator", "reporter_ctry", "partner_ctry", "product_code", "year", "tariff"])

    transform!(df, ["indicator", "reporter_ctry", "partner_ctry"] .=> ByRow(string), renamecols=false)
    df.partner_ctry .= ifelse.(df.indicator .== "MFN", "WTO_MFN", df.partner_ctry)

    # use MFN tariff for tariffs which are missing
    df_MFN = subset(df, :indicator => x -> x .== "MFN")
    df_BPT = subset(df, :indicator => x -> x .== "BPT")

    df_comb = leftjoin(df_BPT, df_MFN[:,[:product_code, :tariff]], on=:product_code, makeunique=true) # possibly reversing join would yield more observations
    # but at least we would not lose the tariffs for product aggregates
    rename!(df_comb, ["tariff", "tariff_1"] .=> ["BPT", "MFN"])
    df_comb.tariff = ifelse.(ismissing.(df_comb.BPT), df_comb.MFN, df_comb.BPT)

    df_formated = df_comb[:, Not(["indicator", "BPT", "MFN"])]

    non_zero = round((1 - count(df_formated.tariff .== 0.0)/nrow(df_formated))*100, digits=2)
    println(" - $non_zero percent are non-zero tariffs.")


    return df_formated
end

for i in ctry
    path = dir * "WTO/" * i * "_2018_MFN_BPT.csv"
    df = transform_wto_data(path)
end




