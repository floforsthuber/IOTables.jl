# Script to import data from RData format ".RData"


using DataFrames, RData


years = 2000:2014 # vector with years covered by WIOD Rev. 2016
N = 44 # number of countries 
S = 56 # number of industries
dir = "C:/Users/u0148308/Desktop/raw/" # location of raw data

# -------------- Import WIOD IO Tables from RData file ----------------------------------------------------------------------------------------------------------------


"""
    import_R(dir::String, year::Integer)

The function loads the WIOD Input-Output Table of specified year as DataFrame into the environment.

# Arguments
- `dir::String`: directory of raw data.
- `year::Integer`: specifies the year of the WIOD table which should be imported.

# Output
- `df::DataFrame{Float64}`: WIOT table as DataFrame in wide format.

# Examples
```julia-repl
julia> import_R(dir, 2012)
```
"""
function import_R(dir::String, year::Integer)

    # load data for specified year and transform into a DataFrame (automatically)
    path = dir * "WIOT" * string(year) * "_October16_ROW.RData"
    df = RData.load(path)["wiot"] # 2472×2690 DataFrame
    transform!(df, [:Year, :RNr] .=> ByRow(Int64) .=> [:Year, :RNr], renamecols=false) # Years, row-industry-identifier as type: Int64

    return df
end



"""
    create_Z_VA(df::DataFrame, N::Integer, S::Integer)

The function creates the origin country-industry destination country-industry intermediate demand matrix Z and 
    the country-industry VA vector.

# Arguments
- `df::DataFrame`: specifies the DataFrame used for creating Z.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `Z::Matrix{Float64}`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `VA::Vector{Float64}`: NS×1, country-industry value added vector VA.

# Examples
```julia-repl
julia> create_Z_VA(name_of_dataframe, N, S)
```
"""
function create_Z_VA(df::DataFrame, N::Integer, S::Integer)
    # Rows and columns are already alphabetically-numerically sorted according to country-industry
    # only need to select the right matrix:
        # do not take first 5 columns: country, sector, year description
        # do not take any final demand or output columns
        # do not take any summary rows (last 8 rows)
    Z = df[1:N*S, 6:N*S+5] 
    Z = Matrix(convert.(Float64, Z)) # NS×NS
    Z = ifelse.(Z .< 0.0, 0.0, Z) # take out negative values
    Z = ifelse.(isnan.(Z), 0.0, Z) # take out NaN

    VA = df[2470, 6:N*S+5] # NS×1
    VA = collect(VA) # method to transform DataFrameRows!
    VA = ifelse.(VA .< 0.0, 0.0, VA) # take out negative values
    VA = ifelse.(isnan.(VA), 0.0, VA) # take out NaN

    return Z, VA
end



"""
    create_F_IV(df::DataFrame, N::Integer, S::Integer)

The function creates the origin country-industry destination country final demand matrix F and
extracts the origin country-industry destination country inventory change vector IV used for inventory adjustments.

# Arguments
- `df::DataFrame`: specifies the DataFrame used for creating F.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `F::Matrix{Float64}`: NS×N, origin country-industry destination country final demand matrix F.
- `IV::Matrix{Float64}`: NS×N, origin country-industry destination country inventory change matrix IV.

# Examples
```julia-repl
julia> create_F_IV(name_of_dataframe, N, S)
```
"""
function create_F_IV(df::DataFrame, N::Integer, S::Integer)
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
    inventory_adjustment(Z::Matrix, F::Matrix, IV::Vector, N::Integer, S::Integer)

The function performs the coutry-industry inventory adjustments.

# Arguments
- `Z::Matrix`: origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: origin country-industry destination country final demand matrix F.
- `IV::Vector`: origin country-industry inventory vector IV.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `Z_adj::Matrix{Float64}`: NS×NS, inventory adjusted origin country-industry destination country-industry intermediate demand matrix Z.
- `F_adj::Matrix{Float64}`: NS×N, inventory adjusted origin country-industry destination country final demand matrix F.
- `Y_adj::Vector{Float64}`: NS×1, inventory adjusted country-industry gross output vector Y.

# Examples
```julia-repl
julia> inventory_adjustment(Z, F, IV, N, S)
```
"""
function inventory_adjustment(Z::Matrix, F::Matrix, IV::Matrix, N::Integer, S::Integer)

    # destination country-industry gross output vector
    Y = [sum(Z[i,:]) + sum(F[i,:] + IV[i,:]) for i in 1:N*S] # NS×1
    Y = ifelse.(abs.(Y) .< 1e-6, 0.0, Y) # Antras and Chor (2018) use the same adjustment

    # destination country-industry IV adjustment matrix
    IV_adj = [Y[i]/(Y[i]-IV[i,j]) for i in 1:N*S, j in 1:N] # NS×N
    IV_adj = ifelse.(isnan.(IV_adj), 1.0, IV_adj) # some entries of Y and IV are equal to zero => NaN
    IV_adj = ifelse.(IV_adj .< 0.95, 0.95, IV_adj) # upper/lower bound on inventory changes (i.e. inventory change can only be ±5%)
    IV_adj = ifelse.(IV_adj .> 1.05, 1.05, IV_adj)

    Z_adj = [Z[i, j]*IV_adj[i, ceil(Int, j/S)] for i in 1:N*S, j in 1:N*S] # N*S×N*S, need to spread IV adjustment of country to country-sector
    F_adj = [F[i, j]*IV_adj[i, j] for i in 1:N*S, j in 1:N] # N*S×N
    Y_adj = [sum(Z_adj[i,:]) + sum(F_adj[i,:]) for i in 1:N*S] # NS×1

    # Y does not exactly match Y_adj => could not find out why this is the case

    return Z_adj, F_adj, Y_adj
end



"""
    create_expenditure_shares(Z::Matrix, F::Matrix, Y::Vector, N::Integer, S::Integer)

The function computes the country-industry VA vector respecting inventory adjustments.

# Arguments
- `Z::Matrix`: origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: origin country-industry destination country final demand matrix F.
- `Y::Vector`: origin country-industry gross output vector Y.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `γ::Matrix{Float64}`: S×N, country-industry intermediate input expenditure share matrix γ.
- `α::Matrix{Float64}`: S×N, country-industry final expenditure expenditure share matrix α.
- `VA_coeff::Vector{Float64}`: NS×1, country value added coefficients vector VA_coeff.
- `TB_ctry::Vector{Float64}`: NS×1, country trade balance (X-M) vector TB.
- `VA_ctry::Vector{Float64}`: N×1, country value added vector VA_ctry.

# Examples
```julia-repl
julia> create_expenditure_shares(Z, F, Y, VA, N, S)
```
"""
function create_expenditure_shares(Z::Matrix, F::Matrix, Y::Vector, N::Integer, S::Integer)

    # reporting industry destination country-industry level input expenditure shares
    # sum over reporting country only (eq. 37) => sum over column gives total intermediate consumption per country-industry (+ VA = Y)
    sZ = zeros(S, N*S) # S×NS
    for i in 1:S
        for j in 1:N*S
            for k in 0:S:N*S-1
            sZ[i, j] += Z[k+i,j]
            end
        end
    end

    # Country-industry level input expenditure shares
    γ = [sZ[i,j]/Y[j] for i in 1:S, j in 1:N*S] # S×NS, divide by total output (i.e. includes final demand)?
    γ = ifelse.(isnan.(γ), 0.0, γ) # some entries in Y equal to zero => NaN => assume 0

    # Country-industry level value added - VVi in Antras and Chor (2018)
    # [total output - sum of inputs of own sector], i.e. computed VA => respects inventory adjustment
    VA = [Y[j] - sum(Z[:,j]) for j in 1:N*S] # NS×1
    VA = ifelse.(VA .< 0.0, 0.0, VA)
    VA = ifelse.(isnan.(VA), 0.0, VA)

    # Country-industry level value added coefficients
    VA_coeff = VA ./ Y # NS×1
    VA_coeff = ifelse.(isnan.(VA_coeff), 1.0, VA_coeff) # some entries in Y equal to zero => NaN, assume VA coefficient to be 1


    # Country level trade balance (WRONG! => if correct α below is the same with both calculations)
    # used computation in Antras and Chor (2018) quite different
    # need to calculate country level exports/imports (row/columns) (i.e. aggregate and remove intra-country trade)

    E_Z = [sum(Z[i:i+S-1,:]) for i in 1:S:N*S] .- [sum(Z[i:i+S-1,i:i+S-1]) for i in 1:S:N*S] # N×1
    E_F = [sum(F[i:i+S-1,:]) for i in 1:S:N*S] .- [sum(F[i:i+S-1,ceil(Int, i/S)]) for i in 1:S:N*S] # N×1
    E = E_Z .+ E_F # N×1

    M_Z = [sum(Z[:,i:i+S-1]) for i in 1:S:N*S] .- [sum(Z[i:i+S-1,i:i+S-1]) for i in 1:S:N*S] # N×1
    M_F = [sum(F[:,i]) for i in 1:N] .- [sum(F[i,i]) for i in 1:N] # N×1
    M = M_Z .+ M_F # N×1

    TB_ctry = E .- M # N×1

    # Country level value added
    VA_ctry = [sum(VA[i:i+S-1]) for i in 1:S:N*S] # N×1

    # Final demand at country level
    F = ifelse.(iszero.(F), 1e-18, F) # Antras and Chor (2018)
    F_ctry = [sum(F[i:i+S-1,:]) for i in 1:S:N*S] # N×1

    # Country-industry final consumption expenditure shares
    # sum over reporting country
    sF = zeros(S, N) # S×N
    for i in 1:S
        for j in 1:N
            for k in 0:S:N*S-1
                sF[i, j] += F[k+i,j]
            end
        end
    end

    # do not perfectly sum to 1 (sometimes quite big difference 0.66 vs 1?)
    α = [sF[i,j]/F_ctry[j] for i in 1:S, j in 1:N] # S×N, as in code of Antras and Chor (2018)
    #α = [sF[i,j]/(VA_ctry[j]+TB_ctry[j]) for i in 1:S, j in 1:N] # S×N, as in equation (38) of Antras and Chor (2018)
    # results are not the same!

    return γ, α, VA_coeff, TB_ctry, VA_ctry, F_ctry
end

df = import_R(dir, 2014)
Z, VA = create_Z_VA(df, N, S)
F, IV = create_F_IV(df, N, S)
Z, F, Y = inventory_adjustment(Z, F, IV, N, S)

γ, α, VA_coeff, TB_ctry, VA_ctry, F_ctry = create_expenditure_shares(Z, F, Y, N, S)


n = [sum(α[:,j]) for j in 1:N]

m = VA_ctry .- F_ctry
TB_ctry ≈ m