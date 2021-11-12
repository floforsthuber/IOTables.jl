# Script to import data from RData format ".RData"


using DataFrames, RData



years = 2000:2014 # vector with years covered by WIOD Rev. 2016
N = 44 # number of countries 
S = 56 # number of industries


# -------------- Import WIOD IO Tables from RData file ----------------------------------------------------------------------------------------------------------------


"""
    import_R(year::Integer)

The function loads the WIOD Input-Output Table of specified year as DataFrame into the environment.

# Arguments
- `year::Integer`: specifies the year of the WIOD table which should be imported.

# Examples
```julia-repl
julia> import_R(2012)
DF
```
"""
function import_R(year::Integer)

    # load data for specified year and transform into a DataFrame (automatically)
    path = "C:/Users/u0148308/Desktop/raw/WIOT"*string(year)*"_October16_ROW.RData"
    df = RData.load(path)["wiot"] # 2472×2690 DataFrame
    transform!(df, [:Year, :RNr] .=> ByRow(Int64) .=> [:Year, :RNr], renamecols=false) # Years, row-industry-identifier as type: Int64

    return df
end

"""
    create_Z(df::DataFrame)

The function creates the origin country-industry destination country-industry intermediate demand matrix Z.

# Arguments
- `df::DataFrame`: specifies the DataFrame used for creating Z.

# Output
- `Z::Matrix{Float64}`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.

# Examples
```julia-repl
julia> create_Z(name_of_dataframe)

```
"""
function create_Z(df::DataFrame)
    # Rows and columns are already alphabetically-numerically sorted according to country-industry
    # only need to select the right matrix:
        # do not take first 5 columns: country, sector, year description
        # do not take any final demand or output columns
        # do not take any summary rows (last 8 rows)
    Z = df[1:N*S, 6:N*S+5] 
    Z = Matrix(convert.(Float64, Z)) # NS×NS
    Z = ifelse.(Z .< 0.0, 0.0, Z) # take out negative values
    Z = ifelse.(isnan.(Z), 0.0, Z) # take out NaN

    return Z
end

"""
    create_F_IV(df::DataFrame)

The function creates the origin country-industry destination country final demand matrix F and
extracts the origin country-industry destination country inventory change vector IV used for inventory adjustments.

# Arguments
- `df::DataFrame`: specifies the DataFrame used for creating F.

# Output
- `F::Matrix{Float64}`: NS×N, origin country-industry destination country final demand matrix F.
- `IV::Matrix{Float64}`: NS×N, origin country-industry destination country inventory change matrix IV.

# Examples
```julia-repl
julia> create_F_IV(name_of_dataframe)

```
"""
function create_F_IV(df::DataFrame)
    # Rows and columns are already alphabetically-numerically sorted according to country-industry
    # only need to select the right matrix:
        # there are five sources of final demand (households, NGOs, government, capital formation, inventories)
    F = df[1:N*S, N*S+6:end-1] # NS×5*N

    # subset inventory changes for inventory adjustment later on
    inventory_columns = 5:5:5*N
    IV = F[:, inventory_columns]
    IV = Matrix(convert.(Float64, IV)) # NS×N
    IV = ifelse.(isnan.(IV), 0.0, IV) # take out NaN (negative values possible!)
    
    # sum four sources of final demand households, NGOs, government, capital formation
    F = F[:, Not(inventory_columns)] # NS×4*N
    F = Matrix(convert.(Float64, F)) # cannot perform comprehensions on DataFrame
    F = [sum(F[i, j:j+3]) for i in 1:N*S, j in 1:4:4*N] # NS×N
    F = ifelse.(F .< 0.0, 0.0, F) # take out negative values
    F = ifelse.(isnan.(F), 0.0, F) # take out NaN

    return F, IV
end



"""
    IV_adjustment(Z::Matrix, F::Matrix, IV::Vector)

The function performs the coutry-industry inventory adjustments.

# Arguments
- `Z::Matrix`: origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: origin country-industry destination country final demand matrix F.
- `IV::Vector`: origin country-industry ventory vector IV.

# Examples
```julia-repl
julia> IV_adjustment(Z, F, IV)

```
"""
function IV_adjustment(Z::Matrix, F::Matrix, IV::Matrix)

    # destination country-industry gross output vector
    Y = [sum(Z[i,:]) + sum(F[i,:] + IV[i,:]) for i in 1:N*S] # NS×1

    # destination country-industry IV adjustment matrix
    IV_adj = [Y[i]/(Y[i]-IV[i,j]) for i in 1:N*S, j in 1:N] # NS×N

    Z_adj = [Z[i, j]*IV_adj[i, ceil(Int, j/S)] for i in 1:N*S, j in 1:N*S] # N*S×N*S, need to spread IV adjustment of country to country-sector
    F_adj = [F[i, j]*IV_adj[i, j] for i in 1:N*S, j in 1:N] # N*S×N
    Y_adj = [sum(Z_adj[i,:]) + sum(F_adj[i,:]) for i in 1:N*S] # NS×1

    return Y, Y_adj
end




a = import_R(2014)
b = create_Z(a)
c, d = create_F_IV(a)

e, f = IV_adjustment(b, c, d)

g = IV_adjustment(b, c, d)