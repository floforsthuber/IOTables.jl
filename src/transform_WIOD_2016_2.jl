# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script with functions to import and transform raw data
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


# --------------- Import data ---------------------------------------------------------------------------------------------------------------------------------


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
julia> WIOD_2012 = import_R(dir, 2012)
```
"""
function import_R(dir::String, year::Integer)

    # load data for specified year and transform into a DataFrame (automatically)
    path = dir * "WIOT" * string(year) * "_October16_ROW.RData"
    df = RData.load(path)["wiot"] # 2472×2690 DataFrame
    transform!(df, [:Year, :RNr] .=> ByRow(Int64) .=> [:Year, :RNr], renamecols=false) # Years, row-industry-identifier as type: Int64

    return df
end


# --------------- Transformations -----------------------------------------------------------------------------------------------------------------------------


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
julia> Z, VA = create_Z_VA(WIOD_2012, N, S)
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
    VA = collect(VA) # method to transform DataFrameRows into a vector!
    VA = ifelse.(VA .< 0.0, 0.0, VA) # take out negative values
    VA = ifelse.(isnan.(VA), 0.0, VA) # take out NaN

    return Z, VA
end


# ---------------


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
julia> F, IV = create_F_IV(name_of_dataframe, N, S)
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


# ---------------


"""
    inventory_adjustment(Z::Matrix, F::Matrix, IV::Vector, N::Integer, S::Integer)

The function performs the coutry-industry inventory adjustments.

# Arguments
- `Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: NS×N, origin country-industry destination country final demand matrix F.
- `IV::Vector`: NS×N, country-industry destination country inventory change matrix IV.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `Z_adj::Matrix{Float64}`: NS×NS, inventory adjusted origin country-industry destination country-industry intermediate demand matrix Z.
- `F_adj::Matrix{Float64}`: NS×N, inventory adjusted origin country-industry destination country final demand matrix F.
- `Y_adj::Vector{Float64}`: NS×1, inventory adjusted country-industry gross output vector Y.

# Examples
```julia-repl
julia> Z, F, Y = inventory_adjustment(Z, F, IV, N, S)
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

    # adjustments done for zero entries, Antras and Chor (2018)
    Z_adj = ifelse.(iszero.(Z_adj), 1e-18, Z_adj) # NS×NS
    F_adj = ifelse.(iszero.(F_adj), 1e-18, F_adj) # NS×N

    Y_adj = [sum(Z_adj[i,:]) + sum(F_adj[i,:]) for i in 1:N*S] # NS×1, use computed gross output vector
    # Y does not exactly match Y_adj ---- WHY??

    return Z_adj, F_adj, Y_adj
end


# ---------------


"""
    create_trade_shares(Z::Matrix, F::Matrix, Y::Vector, N::Integer, S::Integer)

The function computes the origin country-industry destination country intermediate goods trade share matrix π_Z and 
    origin country-industry destination country final goods trade share matrix π_F.

# Arguments
- `Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: NS×N, origin country-industry destination country final demand matrix F.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `π_Z::Matrix{Float64}`: NS×NS, origin country-industry destination country-industry intermediate goods trade share matrix π_Z.
- `π_F::Matrix{Float64}`: NS×N, origin country-industry destination country foods goods trade share matrix π_Z.

# Examples
```julia-repl
julia> π_Z, π_F = create_trade_shares(Z, F, N, S)
```
"""
function create_trade_shares(Z::Matrix, F::Matrix, N::Integer, S::Integer)

    # Origin country-industry destination country-industry intermediate demand trade (export) shares
    Z_agg = [sum(Z[i,:]) for i in 1:N*S] # NS×1
    π_Z = [Z[i,j]/Z_agg[i] for i in 1:N*S, j in 1:N*S] # NS×NS, sum over columns = 1
    π_Z = ifelse.(isnan.(π_Z), 0.0, π_Z) # should never be the case (due to assumptions above)

    # Origin country-industry destination country final demand trade (export) shares
    F_agg = [sum(F[i,:]) for i in 1:N*S] # NS×1
    π_F = [F[i,j]/F_agg[i] for i in 1:N*S, j in 1:N] # NS×N, sum over columns = 1
    π_F = ifelse.(isnan.(π_F), 0.0, π_F) # should never be needed due to assumptions in previous functions

    ## Antras and Chor (2018)
    Z_agg = [sum(Z[s:S:(N-1)*S+s,j]) for s in 1:S, j in 1:N*S] # S×NS
    π_Z = Z ./ repeat(Z_agg, N) # NS×NS
    π_Z = ifelse.(isnan.(π_Z), 0.0, π_Z)
    
    F_agg = [sum(F[s:S:(N-1)*S+s,j]) for s in 1:S, j in 1:N] # S×N
    π_F = F ./ repeat(F_agg, N) # NS×N
    π_F = ifelse.(isnan.(π_F), 0.0, π_F)


    return π_Z, π_F
end


# ---------------


"""
    create_expenditure_shares(Z::Matrix, F::Matrix, Y::Vector, N::Integer, S::Integer)

The function computes the country value added vector VA_ctry respecting inventory adjustments, intermediate goods expenditure share matrix γ 
    and final goods import expenditure share matrix α. Additionally, country-industry value added coeffiecient matrix VA_coeff as well as 
    country trade balance vector TB_ctry and country final goods import vector F_ctry are produced.

# Arguments
- `Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: NS×N, origin country-industry destination country final demand matrix F.
- `Y::Vector`: NS×1, origin country-industry gross output vector Y.
- `π_Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate goods trade share matrix π_Z.
- `π_F::Matrix`: NS×N, origin country-industry destination country foods goods trade share matrix π_Z.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `γ::Matrix{Float64}`: S×NS, country-industry intermediate input expenditure share matrix γ.
- `α::Matrix{Float64}`: S×N, country-industry final expenditure expenditure share matrix α.
- `VA_coeff::Vector{Float64}`: NS×1, country-industry value added coefficients vector VA_coeff.
- `TB_ctry::Vector{Float64}`: NS×1, country trade balance (X-M) vector TB.
- `VA_ctry::Vector{Float64}`: N×1, country value added vector VA_ctry.
- `F_ctry::Vector{Float64}`: N×1, country final import demand vector VA_ctry.

# Examples
```julia-repl
julia> γ, α, VA_coeff, TB_ctry, VA_ctry, F_ctry = create_expenditure_shares(Z, F, Y, π_Z, π_F, VA, N, S)
```
"""
function create_expenditure_shares(Z::Matrix, F::Matrix, Y::Vector, π_Z::Matrix, π_F::Matrix, N::Integer, S::Integer)

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
    γ = [sZ[i,j]/Y[j] for i in 1:S, j in 1:N*S] # S×NS
    γ = ifelse.(isnan.(γ), 0.0, γ) # some entries in Y equal to zero => NaN => assume 0

    # Country-industry level value added - VVi in Antras and Chor (2018)
    # [total output - sum of inputs of own sector], i.e. computed VA => respects inventory adjustment
    VA = [Y[j] - sum(Z[:,j]) for j in 1:N*S] # NS×1
    VA = ifelse.(VA .< 0.0, 0.0, VA)
    VA = ifelse.(isnan.(VA), 0.0, VA)

    # Country-industry level value added coefficients
    VA_coeff = VA ./ Y # NS×1
    VA_coeff = ifelse.(isnan.(VA_coeff), 1.0, VA_coeff) # some entries in Y equal to zero => NaN, assume VA coefficient to be 1

    # Country level value added
    VA_ctry = [sum(VA[i:i+S-1]) for i in 1:S:N*S] # N×1

    # Final import demand at country level (sum over rows not columns - would give export demand)
    F_ctry = [sum(F[:,j]) for j in 1:N] # N×1

    # Country-industry final consumption expenditure (import) shares
    # sum over reporting country
    sF = zeros(S, N) # S×N
    for i in 1:S
        for j in 1:N
            for k in 0:S:N*S-1
                sF[i, j] += F[k+i,j]
            end
        end
    end

    # Computation as in code of Antras and Chor (2018), columns sum to 1
    α = [sF[i,j]/F_ctry[j] for i in 1:S, j in 1:N] # S×N
    #α = [sum(F[j:S:(i-1)*S+j,i])/F_ctry[i] for i in 1:N, j in 1:S] # S×N, Antras and Chor (2018)
    #α2 = [sF[i,j]/(VA_ctry[j] + TB_ctry[j]) for i in 1:S, j in 1:N] # S×N, computation as in equation (38) of Antras and Chor (2018)
    # although F_ctry .== VA_ctry .+ TB_ctry holds true for most countries,
    # α2 considerably different (does not sum to 1) ---- WHY? (even for n=4 for which difference=0)

    # Country level trade balance (exports - imports)
    # used computation in Antras and Chor (2018) quite different
    # need to calculate country level exports/imports (row/columns) (i.e. aggregate and remove intra-country trade)
    E = [sum(Y[i:i+S-1]) - sum(Z[i:i+S-1,i:i+S-1]) - sum(F[i:i+S-1,ceil(Int, i/S)]) for i in 1:S:N*S] # N×1

    M_Z = [sum(Z[:,j:j+S-1]) - sum(Z[j:j+S-1,j:j+S-1]) for j in 1:S:N*S]
    M_F = [sum(F[:,ceil(Int, j/S)]) - sum(F[j:j+S-1,ceil(Int, j/S)]) for j in 1:S:N*S]
    M = M_Z .+ M_F # N×1

    #### Antras and Chor (2018) calculation
    E_A = zeros(N)
    M_A = zeros(N)
    for j in 1:N
        for s in 1:N
            for i in 1:N
                iiss = (i-1)*S+s
                jjss = (j-1)*S+s
                for r in 1:S
                    jjrr = (j-1)*S+r
                    iirr = (i-1)*S+r
                    M_A[j] += π_Z[iiss,jjrr]*γ[s,jjrr]*Y[jjrr]
                    E_A[j] += π_Z[jjss,iirr]*γ[s,iirr]*Y[iirr]
                end
                M_A[j] += π_F[iiss,j]*α[s,j]*F_ctry[j]
                E_A[j] += π_F[jjss,i]*α[s,i]*F_ctry[i]
            end
        end
    end

    TB_ctry = E .- M # N×1
    TB_ctry = E_A .- M_A # N×1

    return γ, α, VA_coeff, TB_ctry, VA_ctry, F_ctry
end


# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Summary function
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
    transform_WIOD_2016(dir::String, year::Integer)

The function imports and transforms the raw data for further use in the model.

# Arguments
- `dir::String`: directory of raw data.
- `year::Integer`: specifies the year of the WIOD table which should be imported.

# Output
- `Z::Matrix{Float64}`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix{Float64}`: NS×N, origin country-industry destination country final demand matrix F.
- `Y::Vector{Float64}`: NS×1, origin country-industry gross output vector Y.
- `F_ctry::Vector{Float64}`: N×1, country final import demand vector VA_ctry.
- `TB_ctry::Vector{Float64}`: NS×1, country trade balance (X-M) vector TB.
- `VA_ctry::Vector{Float64}`: N×1, country value added vector VA_ctry.
- `VA_coeff::Vector{Float64}`: NS×1, country-industry value added coefficients vector VA_coeff.
- `γ::Matrix{Float64}`: S×NS, country-industry intermediate input expenditure share matrix γ.
- `α::Matrix{Float64}`: S×N, country-industry final expenditure expenditure share matrix α.
- `π_Z::Matrix{Float64}`: NS×NS, origin country-industry destination country-industry intermediate goods trade share matrix π_Z.
- `π_F::Matrix{Float64}`: NS×N, origin country-industry destination country final goods trade share matrix π_F.

# Examples
```julia-repl
julia> Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_WIOD_2016(dir, 2014)
```
"""
function transform_WIOD_2016(dir::String, year::Integer)

    df = import_R(dir, year)
    Z, VA = create_Z_VA(df, N, S)
    F, IV = create_F_IV(df, N, S)
    Z, F, Y = inventory_adjustment(Z, F, IV, N, S)
    π_Z, π_F = create_trade_shares(Z, F, N, S)
    γ, α, VA_coeff, TB_ctry, VA_ctry, F_ctry = create_expenditure_shares(Z, F, Y, π_Z, π_F, N, S)


    return Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F
end