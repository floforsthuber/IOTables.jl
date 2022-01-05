# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Baseline model
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData,  XLSX, LinearAlgebra, Statistics, CSV

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

df_tariffs = DataFrame(reporter_lab=String[], year=Int64[], HS2017_code=String[], 
    HS2017_lab=String[], duty_lab=String[], avg_tariff=Union{Float64,Missing}[], reporter=String[])

function tariff_data(reporter::String, id::String, year::String)

    col_names_old = ["Reporter", "Year", "Product Code", "Product Description", "Duty Description", "All items excl. NA - Simple Averages"]
    col_names_new = ["reporter_lab", "year", "HS2017_code", "HS2017_lab", "duty_lab", "avg_tariff"]

    path = dir * "MFN/" * reporter * "_" * year * "_TA/" * reporter * "_" * year * "_Duties_TA.txt"

    df = CSV.read(path, DataFrame)
    df = df[:, col_names_old]
    rename!(df, col_names_new)
    df[:, :reporter] .= id

    return df
end

ctry_names_2013 = ["AUS", "AUT", "BEL", "BGR", "BRA", "CAN", "CHN", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "GBR", 
    "GRC", "HUN", "IDN", "IND", "IRL", "ITA", "JPN", "KOR", "LTU", "LUX", "LVA", "MEX", "MLT", "NLD", "POL", "PRT", "ROM", "RUS", 
    "SVK", "SVN", "SWE", "TUR", "TWN", "USA"]

ctry_EU28 = ["AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GBR", "GRC", "HUN", "IRL", "ITA", "LVA", 
    "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROM", "SVK", "SVN", "ESP", "SWE"]


for i in ctry_names_2013

    if i in ctry_EU28 
        df = tariff_data("EU", i, "2018") # for EU countries
    elseif i in ["RUS", "TUR"] 
        df = tariff_data(i, i, "2019") # countries for which tariffs in 2018 were not available
    else
        df = tariff_data(i, i, "2018")
    end

    append!(df_tariffs, df)
end


# Correspondence table HS2017, ISIC rev. 3 and ISIC rev. 4
correspondence_HS_ISIC = DataFrame(XLSX.readtable(dir * "correspondence_HS_ISIC.xlsx", "Sheet1")...)
transform!(correspondence_HS_ISIC, [:ISIC3, :ISIC4, :HS2017] .=> ByRow(x-> lpad(x,2,'0')), renamecols=false)
rename!(correspondence_HS_ISIC, [:HS2017, :ISIC3, :ISIC4] .=> [:HS2017_code, :ISIC3_code, :ISIC4_code])

# join tables
transform!(df_tariffs, :HS2017_code => ByRow(x->replace(x, "'" => "")), renamecols=false) # remove ' from column otherwise we cannot join properly
df = leftjoin(df_tariffs, correspondence_HS_ISIC, on=:HS2017_code)

# -------


# WIOD rev. 2013
gdf = groupby(df, [:reporter, :ISIC3_code]) # group to compute average tariffs according to ISIC rev. 3 classification
df_tariffs = combine(gdf, :avg_tariff => mean, renamecols=false)
sort!(df_tariffs, [:reporter, :ISIC3_code])
dropmissing!(df_tariffs) # drop missing for using average later

# use average of all countries for RoW
gdf = groupby(df_tariffs, :ISIC3_code)
df_RoW = combine(gdf, :avg_tariff => mean, renamecols=false)
df_RoW[:, :reporter] .= "ZoW"

df_all = vcat(df_tariffs, df_RoW)



reporter = [ctry_names_2013; "ZoW"] # use "ZoW" instead of "RoW" so we can sort later on
industries = lpad.(1:S, 2, '0') # left pad with leading zeros so we can sort later on

gdf = groupby(df_all, :reporter)
df_ctry_avg = combine(gdf, :avg_tariff => mean, renamecols=false)
avg_ctry_tariff = df_ctry_avg.avg_tariff # use average country tariffs for sectors which are missing

df = DataFrame([repeat(reporter, inner=S) repeat(industries, outer=N) repeat(avg_ctry_tariff, inner=S)], [:reporter, :ISIC3_code, :ctry_avg])
df = leftjoin(df, df_all, on=[:reporter, :ISIC3_code])
df.avg_tariff = ifelse.(ismissing.(df.avg_tariff) .== true, df.ctry_avg, df.avg_tariff)
rename!(df, :ISIC3_code => :industry)
df_tariffs = df[:, [:reporter, :industry, :avg_tariff]]
sort!(df_tariffs, [:reporter, :industry])

CSV.write(dir * "MFN/df_tariffs.csv", df_tariffs)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# format data to NS×N matrix
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

vec_tariffs = df_tariffs.avg_tariff
matrix_tariffs = repeat(vec_tariffs, 1, N)
countries = [ctry_names_2013; "ZoW"]

τ_Z = ones(N*S,N)

for (i, ctry_i) in enumerate(countries)
    for (j, ctry_j) in enumerate(countries)

        if ctry_i == ctry_j
            τ_Z[1+(i-1)*S:i*S,j] .= 0

        elseif ctry_i in ctry_EU28 && ctry_j in ctry_EU28
            τ_Z[1+(i-1)*S:i*S,j] .= 0

        else
            τ_Z[1+(i-1)*S:i*S,j] = matrix_tariffs[1+(i-1)*S:i*S,j]
        end

    end
end

# export
rows = repeat(countries, inner=S) .* "__" .* repeat(industries, outer=N)
df_tariff_matrix = DataFrame([rows τ_Z], ["reporter"; countries])
XLSX.writetable(dir * "MFN/tariff_matrix.xlsx", df_tariff_matrix, overwrite=true)