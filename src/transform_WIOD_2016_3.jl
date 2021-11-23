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
    create_matrices(df::DataFrame, N::Integer, S::Integer)

The function creates the origin country-industry destination country-industry intermediate demand matrix Z and 
    the country-industry VA vector.

# Arguments
- `df::DataFrame`: specifies the DataFrame used for creating Z.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `Z::Matrix{Float64}`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix{Float64}`: NS×N, origin country-industry destination country final demand matrix F.
- `Y::Vector{Float64}`: NS×1, inventory adjusted country-industry gross output vector Y.
- `VA::Vector{Float64}`: NS×1, country-industry value added vector VA.

# Examples
```julia-repl
julia> Z, F, Y, VA = create_matrices(WIOD_2012, N, S)
```
"""
function create_matrices(df::DataFrame, N::Integer, S::Integer)

    # Intermediate demand
    Z = df[1:N*S, 6:N*S+5] # NS×NS
    Z = Matrix(convert.(Float64, Z))

    # Final demand
    F = df[1:N*S, N*S+6:end-1] # NS×5*N
    
    # Subset final demand for inventory adjustments
    inventory_columns = 5:5:5*N
    IV = F[:, inventory_columns] # NS×N
    F = F[:, Not(inventory_columns)] # NS×4*N
    
    IV = Matrix(convert.(Float64, IV)) # NS×N
    F = Matrix(convert.(Float64, F)) # NS×4*N

    IV = [sum(IV[i,:]) for i in 1:N*S] # NS×1

    # Compute gross output from intermediate and final demand
    Y = [sum(Z[i,:]) + sum(F[i,:]) + IV[i] for i in 1:N*S] # NS×1
    Y = ifelse.(abs.(Y) .< 1e-6, 0.0, Y)
   
    # Inventory adjustment (spread IV regardless of origin over all columns)
    adj = Y .- IV # NS×1
    adj = ifelse.(adj .< 0.0, 0.0, adj) # in case negative

    Z = Z .* (repeat(Y, 1, N*S) ./ repeat(adj, 1, N*S)) # NS×NS
    Z = ifelse.(isnan.(Z), 0.0, Z)
    Z = ifelse.(isinf.(Z), 0.0, Z)

    F = F .* (repeat(Y, 1, N*4) ./ repeat(adj, 1, N*4)) # NS×N*4, inventory adjustment
    F = ifelse.(isnan.(F), 0.0, F)
    F = ifelse.(isinf.(F), 0.0, F)
    F = [sum(F[i,j:j+3]) for i in 1:N*S, j in 1:4:N*4] # NS×N

    # Use computed values to have internal consistency
    Y = [sum(Z[i,:]) + sum(F[i,:]) for i in 1:N*S] # NS×1
    VA = [Y[i] - sum(Z[:,i]) for i in 1:N*S] # NS×1, compute value added (gross output minus inputs)

    return Z, F, Y, VA
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
- `VA::Vector`: NS×1, country-industry value added vector VA.
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
function create_expenditure_shares(Z::Matrix, F::Matrix, Y::Vector, VA::Vector, π_Z::Matrix, π_F::Matrix, N::Integer, S::Integer)

    # Country-industry level input expenditure shares
    # sum over reporting country only (eq. 37) => sum over column gives total intermediate consumption per country-industry (+ VA = Y)
    Z_agg = [sum(Z[i:S:(N-1)*S+i, j]) for i in 1:S, j in 1:N*S]
    Z_agg_VA = [Z_agg; VA']
    γ_VA = [Z_agg_VA[i,j]/sum(Z_agg_VA[:,j]) for i in 1:S+1, j in 1:N*S]
    γ = γ_VA[1:S, 1:N*S] # S×NS
    γ = ifelse.(isnan.(γ), 0.0, γ) # some entries in potentially equal to zero => NaN => assume 0

    # Country-industry level value added coefficients
    VA_coeff = γ_VA[S+1, 1:N*S] # NS×1
    VA_coeff = ifelse.(isnan.(VA_coeff), 1.0, VA_coeff) # some entries in potentially equal to zero => NaN, assume VA coefficient to be 1

    # Country level value added
    VA_ctry = [sum(VA[i:i+S-1]) for i in 1:S:N*S] # N×1

    # Final import demand at country level (sum over rows not columns - would give export demand)
    F_ctry = [sum(F[:,j]) for j in 1:N] # N×1

    # Country-industry final consumption expenditure (import) shares => columns sum to 1
    α = [sum(F[i:S:(N-1)*S+i, j])/F_ctry[j] for i in 1:S, j in 1:N] # S×N

    # Country level trade balance (exports - imports)
    # used computation in Antras and Chor (2018) quite different
    # need to calculate country level exports/imports (row/columns) (i.e. aggregate and remove intra-country trade)
    E = [sum(Y[i:i+S-1]) - sum(Z[i:i+S-1,i:i+S-1]) - sum(F[i:i+S-1,ceil(Int, i/S)]) for i in 1:S:N*S] # N×1

    M_Z = [sum(Z[:,j:j+S-1]) - sum(Z[j:j+S-1,j:j+S-1]) for j in 1:S:N*S]
    M_F = [sum(F[:,ceil(Int, j/S)]) - sum(F[j:j+S-1,ceil(Int, j/S)]) for j in 1:S:N*S]
    M = M_Z .+ M_F # N×1

    # #### Antras and Chor (2018) calculation
    # E_A = zeros(N)
    # M_A = zeros(N)
    # for j in 1:N
    #     for s in 1:N
    #         for i in 1:N
    #             iiss = (i-1)*S+s
    #             jjss = (j-1)*S+s
    #             for r in 1:S
    #                 jjrr = (j-1)*S+r
    #                 iirr = (i-1)*S+r
    #                 M_A[j] += π_Z[iiss,jjrr]*γ[s,jjrr]*Y[jjrr]
    #                 E_A[j] += π_Z[jjss,iirr]*γ[s,iirr]*Y[iirr]
    #             end
    #             M_A[j] += π_F[iiss,j]*α[s,j]*F_ctry[j]
    #             E_A[j] += π_F[jjss,i]*α[s,i]*F_ctry[i]
    #         end
    #     end
    # end

    TB_ctry = E .- M # N×1
    #TB_ctry = E_A .- M_A # N×1

    #TB_ctry == VA_ctry .- F_ctry # must hold, holds much better with own calculation?

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
function transform_WIOD_2016(dir::String, year::Integer, N::Integer, S::Integer)

    df = import_R(dir, year)
    Z, F, Y, VA = create_matrices(df, N, S)
    π_Z, π_F = create_trade_shares(Z, F, N, S)
    γ, α, VA_coeff, TB_ctry, VA_ctry, F_ctry = create_expenditure_shares(Z, F, Y, VA, π_Z, π_F, N, S)

    return Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F
end