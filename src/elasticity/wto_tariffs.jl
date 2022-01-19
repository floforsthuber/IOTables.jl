# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script to prepare WTO tariff data from 2018
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData, XLSX, LinearAlgebra, Statistics, CSV

# dir = "X:/VIVES/1-Personal/Florian/git/IOTables/src/"
# include(dir * "model/transform_data.jl") # Script with functions to import and transform raw data
# include(dir * "model/price_hat.jl") # Script with function to obtain the price index
# include(dir * "model/wage_hat.jl") # Script with function to obtain the wages and gross output
# include(dir * "counterfactual/tariffs_function.jl") # Script with functions to create τ_hat_Z, τ_hat_F from tariff data
# include(dir * "counterfactual/head_ries_index.jl") # Script with functions to create bilateral Head-Ries index (symmetric bilateral trade costs)
# include(dir * "elasticity/elasticities_functions.jl") # Script with functions to compute statistics for estimating sectoral trade elasticities (method of Caliendo and Parro, 2015)

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

year = 2018 # year for tariff data

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# prepare tariff data from WTO for EU
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

path = dir * "WTO/EU.csv"

df2 = CSV.read(path, DataFrame, delim=',', quoted=true, quotechar='"')

# problem with reading data (one column gets split into multiple)
cols = ["Indicator", "Reporting Economy ISO3A Code", "Partner Economy ISO3A Code", "Product/Sector Code", "Year", "Value", "Column25", "Column26", "Column27", "Column28"]
df = df2[:, cols]

transform!(df, :Column28 => ByRow(x -> ifelse(ismissing(x), missing, string(x))), renamecols=false) # tariffs imported as Int
# remove all strings which are split over columns and fill with NaN instead
transform!(df, [:Value, :Column25, :Column26, :Column27, :Column28] .=> ByRow(x -> ifelse(ismissing(x) || !contains(x, " "), x, NaN)), renamecols=false)

# function to extract the tariff value or missing from the split columns
function extract_tariff_value(a, b, c, d, e)
    v = [a, b, c, d, e] # split columns which contain the tariff values
    index = findall(typeof.(v) .== String)
    value = ifelse(any(typeof.(v) .== String), parse.(Float64, v[index]), [missing])
    return value[1]
end

df.tariff = extract_tariff_value.(df.Value, df.Column25, df.Column26, df.Column27, df.Column28)

# format table
df = df[:, Not(["Value", "Column25", "Column26", "Column27", "Column28"])]
df.Indicator .= ifelse.(df.Indicator .== "HS MFN - Simple average with ad valorem equivalents (AVE)", "MFN", "BPT")
rename!(df, ["indicator", "reporter_ctry", "partner_ctry", "product_code", "year", "tariff"])
transform!(df, ["indicator", "reporter_ctry", "partner_ctry"] .=> ByRow(string), renamecols=false)
df.partner_ctry .= ifelse.(df.indicator .== "MFN", "WTO_MFN", df.partner_ctry)

# subset data for specified/closest year
println(" EU:")
years_with_MFN = unique(subset(df, :indicator => x -> x .== "MFN").year)
years_with_BPT = unique(subset(df, :indicator => x -> x .== "BPT").year)


if (year in years_with_MFN) & (year in years_with_BPT)
    closest_year_MFN = year
    closest_year_BPT = year
    println(" - Tariff data for $year is available for both MFN and BPT.")
else
    if length(years_with_MFN[years_with_MFN .< year]) != 0
        prev_years = years_with_MFN[years_with_MFN .< year] # take data for previous years since those are in effect
        closest_year_MFN = maximum(prev_years) # finds closest year for which data is available
    else
        prev_years = years_with_MFN[years_with_MFN .> year]
        closest_year_MFN = minimum(prev_years)
    end

    if length(years_with_BPT[years_with_BPT .< year]) != 0
        prev_years = years_with_BPT[years_with_BPT .< year] 
        closest_year_BPT = maximum(prev_years)
    else
        prev_years = years_with_BPT[years_with_BPT .> year] 
        closest_year_BPT = minimum(prev_years)
    end

    println(" - Tariff data for $year is not available, taking $closest_year_MFN for MFN and $closest_year_BPT for BPT.")
end

# use MFN tariff for tariffs which are missing
df_MFN = subset(df, :indicator => x -> x .== "MFN", :year => x -> x .== closest_year_MFN)
df_BPT = subset(df, :indicator => x -> x .== "BPT", :year => x -> x .== closest_year_BPT)

# ------

df_t = copy(df_MFN)
ctry = ["AUS", "BRA", "CAN", "CHN", "IDN", "IND", "JPN", "KOR", "RUS", "TUR", "TWN", "USA"] # countries with WTO data
ctry = [ctry; "EEC"] # add RoW and EU

for i in ctry
    df_loop = copy(df_MFN)
    df_loop.partner_ctry .= i
    df_loop.indicator .= "BPT"
    append!(df_t, df_loop)
end

# take BPT tariff if available
df_join = leftjoin(df_t, df_BPT[:, [:partner_ctry, :product_code, :tariff]], on=[:partner_ctry, :product_code], makeunique=true)
rename!(df_join, ["tariff", "tariff_1"] .=> ["MFN", "BPT"])
df_join.tariff = ifelse.(ismissing.(df_join.BPT), df_join.MFN, df_join.BPT)
df_tariffs = df_join[:, Not(["indicator", "BPT", "MFN"])] # new collection of tariffs for countries which use MFN + BPT


# ------

# conversion from HS2017 6digits to ISIC rev. 4 (using HS2012 for now)
df_corr = DataFrame(XLSX.readtable(dir * "HS_ISIC_BEC_EUC.xlsx", "Sheet1")...)

df_corr = df_corr[:, ["HS_6digit", "ISIC_first_2_digits"]]
rename!(df_corr, ["product_code", "ISIC"])
transform!(df_corr, [:product_code, :ISIC] .=> ByRow(string), renamecols=false)

# ------

function aggregation(df::DataFrame, df_corr::DataFrame)
    df.product_code = lpad.(df.product_code, 6, '0')

    df_join = leftjoin(df, df_corr, on=:product_code)
    subset!(df_join, :ISIC => ByRow(x -> !ismissing(x) && x != "missing")) # take out missing values (MFN also includes 2-4-digit codes)
    gdf = groupby(df_join, [:reporter_ctry, :partner_ctry, :ISIC])
    df = combine(gdf, :tariff => (x -> mean(skipmissing(x))) => :tariff)
    sort!(df, [:partner_ctry, :ISIC])

    df.tariff .= ifelse.(df.reporter_ctry .== df.partner_ctry, 0.0, df.tariff)
    df.partner_ctry .= ifelse.(df.partner_ctry .== "WTO_MFN", "ZoW", df.partner_ctry) # take MFN tariffs for RoW

    return df
end

# ------

df = aggregation(df_tariffs, df_corr)

missing_values = round(count(ismissing.(df.tariff))/nrow(df)*100, digits=2)
non_zero = round((1 - count(skipmissing(df.tariff .== 0.0))/nrow(df))*100, digits=2)
println(" - $missing_values percent are missing. Out of non-missing values, $non_zero percent are non-zero tariffs.")


# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# prepare tariff data from WTO for other countries
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

function transform_wto_data(path::String)
    
    df = CSV.read(path, DataFrame, delim=',', quoted=true, quotechar='"')
    cols = ["Indicator", "Reporting Economy ISO3A Code", "Partner Economy ISO3A Code", "Product/Sector Code", "Year", "Value"]
    df = df[:, cols]

    df.Indicator .= ifelse.(df.Indicator .== "HS MFN - Simple average with ad valorem equivalents (AVE)", "MFN", "BPT")
    rename!(df, ["indicator", "reporter_ctry", "partner_ctry", "product_code", "year", "tariff"])

    transform!(df, ["indicator", "reporter_ctry", "partner_ctry"] .=> ByRow(string), renamecols=false)
    df.partner_ctry .= ifelse.(df.indicator .== "MFN", "WTO_MFN", df.partner_ctry)

    # subset data for specified year

    if length(unique(df.indicator)) == 2 # both MFN and BPT are available

        years_with_MFN = unique(subset(df, :indicator => x -> x .== "MFN").year)
        years_with_BPT = unique(subset(df, :indicator => x -> x .== "BPT").year)

        if (year in years_with_MFN) & (year in years_with_BPT)
            closest_year_MFN = year
            closest_year_BPT = year
            println(" - Tariff data for $year is available for both MFN and BPT.")
        else
            if length(years_with_MFN[years_with_MFN .< year]) != 0
                prev_years = years_with_MFN[years_with_MFN .< year] # take data for previous years since those are in effect
                closest_year_MFN = maximum(prev_years) # finds closest year for which data is available
            else
                prev_years = years_with_MFN[years_with_MFN .> year]
                closest_year_MFN = minimum(prev_years)
            end
        
            if length(years_with_BPT[years_with_BPT .< year]) != 0
                prev_years = years_with_BPT[years_with_BPT .< year] 
                closest_year_BPT = maximum(prev_years)
            else
                prev_years = years_with_BPT[years_with_BPT .> year] 
                closest_year_BPT = minimum(prev_years)
            end
        
            println(" - Tariff data for $year is not available, taking $closest_year_MFN for MFN and $closest_year_BPT for BPT.")
        end

        df_MFN = subset(df, :indicator => x -> x .== "MFN", :year => x -> x .== closest_year_MFN)
        df_BPT = subset(df, :indicator => x -> x .== "BPT", :year => x -> x .== closest_year_BPT)


    else 
        years_with_MFN = unique(subset(df, :indicator => x -> x .== "MFN").year) # at least need to have MFN data

        if year in years_with_MFN
            closest_year_MFN = year
            println(" - Tariff data for $year is available, but only MFN is available.")
        else
            if length(years_with_MFN[years_with_MFN .< year]) != 0
                prev_years = years_with_MFN[years_with_MFN .< year] # take data for previous years since those are in effect
                closest_year_MFN = maximum(prev_years) # finds closest year for which data is available
            else
                prev_years = years_with_MFN[years_with_MFN .> year]
                closest_year_MFN = minimum(prev_years)
            end
    
            println(" - Tariff data for $year is not available, taking $closest_year_MFN for MFN (BPT not available).")
        end

        df_MFN = subset(df, :indicator => x -> x .== "MFN", :year => x -> x .== closest_year_MFN)
        df_BPT = copy(df_MFN) # only have MFN tariffs

    end
    

    # ------

    df_t = copy(df_MFN) # initialize

    # create MFN tariffs for each country in IO Table and then match BPT
    for i in ctry
        df_loop = copy(df_MFN)
        df_loop.partner_ctry .= i
        df_loop.indicator .= "BPT"
        append!(df_t, df_loop)
    end

    # take BPT tariff if available
    df_join = leftjoin(df_t, df_BPT[:, [:partner_ctry, :product_code, :tariff]], on=[:partner_ctry, :product_code], makeunique=true)
    rename!(df_join, ["tariff", "tariff_1"] .=> ["MFN", "BPT"])
    df_join.tariff = ifelse.(ismissing.(df_join.BPT), df_join.MFN, df_join.BPT)
    df_tariffs = df_join[:, Not(["indicator", "BPT", "MFN"])] # new collection of tariffs for countries which use MFN + BPT

    df = aggregation(df_tariffs, df_corr)

    missing_values = round(count(ismissing.(df.tariff))/nrow(df)*100, digits=2)
    non_zero = round((1 - count(skipmissing(df.tariff .== 0.0))/nrow(df))*100, digits=2)
    println(" - $missing_values percent are missing. Out of non-missing values, $non_zero percent are non-zero tariffs.")

    return df
end

wto_ctry = ["AUS", "BRA", "CAN", "CHN", "IDN", "IND", "JPN", "KOR", "MEX", "RUS", "TUR", "TWN", "USA"] # countries in IO table


for i in wto_ctry
    path = dir * "WTO/" * i * ".csv"
    println("\n $i:")
    append!(df, transform_wto_data(path))
end

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# format tariffs to into NS×N matrix
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


df.reporter_ctry .= ifelse.(df.reporter_ctry .== "CHT", "TWN", df.reporter_ctry)
df.partner_ctry .= ifelse.(df.partner_ctry .== "CHT", "TWN", df.partner_ctry)

ctry_names_2013 = ["AUS", "AUT", "BEL", "BGR", "BRA", "CAN", "CHN", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "GBR", 
    "GRC", "HUN", "IDN", "IND", "IRL", "ITA", "JPN", "KOR", "LTU", "LUX", "LVA", "MEX", "MLT", "NLD", "POL", "PRT", "ROM", "RUS", 
    "SVK", "SVN", "SWE", "TUR", "TWN", "USA"]

ctry_EU28 = ["AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GBR", "GRC", "HUN", "IRL", "ITA", "LVA", 
    "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROM", "SVK", "SVN", "ESP", "SWE"]

countries = [ctry_names_2013; "ZoW"]

industries = lpad.(1:S, 2, '0')

τ_Z = ones(N*S,N)

for (i, ctry_i) in enumerate(countries)
    for (s, ind) in enumerate(industries)
        for (j, ctry_j) in enumerate(countries)

            if ctry_i == ctry_j # within country trade
                τ_Z[(i-1)*S+s,j] = 0.0

            elseif ctry_i in ctry_EU28 && ctry_j in ctry_EU28 # within EU trade
                τ_Z[(i-1)*S+s,j] = 0.0

            else
                ctry_i = ifelse(ctry_i in ctry_EU28, "EEC", ctry_i)
                ctry_j = ifelse(ctry_j in ctry_EU28, "EEC", ctry_j)

                # assume zero tariffs for countries/industries which are not available
                # assume RoW imposes zero tariffs on all trading partners
                value = df.tariff[(df.reporter_ctry .== ctry_i) .& (df.partner_ctry .== ctry_j) .& (df.ISIC .== ind)]
                tariff = ifelse(isempty(value), 0.0, value)

                τ_Z[(i-1)*S+s,j] = tariff[1]
            end

        end
    end
end

# export data
rows = repeat(countries, inner=S) .* "__" .* repeat(industries, outer=N)
df_tariff_matrix = DataFrame([rows τ_Z], ["reporter"; countries])
XLSX.writetable(dir * "WTO/tariff_matrix_" * string(year) * ".xlsx", df_tariff_matrix, overwrite=true)