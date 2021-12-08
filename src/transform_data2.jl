# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script with functions to import and transform raw data
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


# --------------- Raw data ---------------------------------------------------------------------------------------------------------------------------------

"""
    import_data(dir::String, revision::String, year::Integer, N::Integer, S::Integer)

The function loads the WIOD Input-Output Table of specified year as DataFrame into the environment and 
    extracts the raw intermediate, final demand matrices (Z, F) and the inventory adjustment matrix IV.

# Arguments
- `dir::String`: directory of raw data (either folder location or entire directory).
- `revision::String`: specifies the WIOD revision (2013 or 2016) as a string.
- `year::Integer`: specifies the year of the WIOD table which should be imported.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `Z::Matrix{Float64}`: NS×NS, raw origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix{Float64}`: NS×N, raw origin country-industry destination country final demand matrix F.
- `IV::Matrix{Float64}`: NS×N, raw origin country-industry destination country inventory adjustment matrix IV.

# Examples
```julia-repl
julia> Z, F, IV = import_data("C:/Users/u0148308/Desktop/raw/", "2016", 2014, 44, 56)
```
"""

function import_data(dir::String, revision::String, year::Integer, N::Integer, S::Integer)

    # switch between WIOD revisions (2016 uses R, 2013 uses Excel as input)
    if revision == "2016"
        file = "WIOT" * string(year) * "_October16_ROW.RData"
        path = ifelse(contains(dir[end-5:end], '.'), dir, dir * file)

        df = RData.load(path)["wiot"]
        df = df[1:N*S, 6:end-1] # NS×NS+5N, take out the rows/columns with country/industry names

    else
        name_2013 = ifelse(year >= 2008, "Sep", "Apr")
        file = "WIOT" * string(year)[end-1:end] * "_ROW_" * name_2013 * "12.xlsx"
        path = ifelse(contains(dir[end-5:end], '.'), dir, dir * file)
    
        df = DataFrames.DataFrame(XLSX.readxlsx(path)["WIOT_$year"][:], :auto)
        df = df[7:end-8,5:end-1] # NS×NS+5N, take out the rows/columns with country/industry names

    end

    IO_table = Matrix(convert.(Float64, df)) # NS×NS+5N
    inventory_columns = N*S+5:5:N*S+5*N
    
    Z = IO_table[:, 1:N*S] # NS×NS
    F = IO_table[:, Not([1:N*S; inventory_columns])] # NS×4N
    IV = IO_table[:, inventory_columns] # NS×N

    println(" ✓ Raw data from WIOD (rev. $revision) for $year was successfully imported!")

    return Z, F, IV
end


# --------------- Transformations -----------------------------------------------------------------------------------------------------------------------------

"""
    inventory_adjustment(Z::Matrix, F::Matrix, IV::Matrix, N::Integer, S::Integer)

The function performs the inventory adjustment procedure proposed in Antras et al. (2012).

# Arguments
- `Z::Matrix`: NS×NS, raw origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: NS×N, raw origin country-industry destination country final demand matrix F.
- `IV::Matrix`: NS×N, raw origin country-industry destination country inventory adjustment matrix IV.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `Z::Matrix{Float64}`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix{Float64}`: NS×N, origin country-industry destination country final demand matrix F.
- `Y::Vector{Float64}`: NS×1, inventory adjusted country-industry gross output vector Y.
- `VA::Vector{Float64}`: NS×1, country-industry value added vector VA.

# Examples
```julia-repl
julia> Z, F, Y, VA = inventory_adjustment(Z, F, IV, N, S)
```
"""

function inventory_adjustment(Z::Matrix, F::Matrix, IV::Matrix, N::Integer, S::Integer)

    # remove zero and negative entries (page 37 Antras and Chor, 2018)
    #Z = ifelse.(Z .<= 0.0, 1e-18, Z) 
    #F = ifelse.(F .<= 0.0, 1e-18, F)

    Z = ifelse.(Z .< 0.0, 0.0, Z) # no negative values possible in Z
    F = ifelse.(F .< 0.0, 0.0, F) # no negative values possible in F

    # Inventory adjustment (spread IV regardless of origin over all columns)
    IV = [sum(IV[i,:]) for i in 1:N*S] # NS×1, negative values possible
    Y = [sum(Z[i,:]) + sum(F[i,:]) + IV[i] for i in 1:N*S] # NS×1, possibly negative if IV is negative enough
    Y = ifelse.(Y .< 0.0, 0.0, Y) # exclude negative values 

    adjustment = Y .- IV # NS×1
    adjustment = ifelse.(adjustment .< 0.0, 0.0, adjustment) # in case negative, as in Antras and Chor (2018)

    # Inventory adjustment
    Z = Z .* ( repeat(Y, 1, N*S) ./ repeat(adjustment, 1, N*S) ) # NS×NS, inventory adjustment
    Z = ifelse.(isnan.(Z), 0.0, Z)
    Z = ifelse.(isinf.(Z), 0.0, Z)
    Z = ifelse.(Z .< 0.0, 0.0, Z) # NS×NS

    F = F .* ( repeat(Y, 1, N*4) ./ repeat(adjustment, 1, N*4) ) # NS×N*4, inventory adjustment
    F = ifelse.(isnan.(F), 0.0, F)
    F = ifelse.(isinf.(F), 0.0, F)
    F = ifelse.(F .< 0.0, 0.0, F)
    F = [sum(F[i,j:j+3]) for i in 1:N*S, j in 1:4:N*4] # NS×N, sum all types of final demand

    # Use computed values for Y and VA to have internal consistency
    Y = [sum(Z[i,:]) + sum(F[i,:]) for i in 1:N*S] # NS×1

    VA = [Y[i] - sum(Z[:,i]) for i in 1:N*S] # NS×1, compute value added (gross output minus sum of inputs)
    #VA = ifelse.(VA .< 0.0, 0.0, VA) # in case one is negative (actually problematic since Z, F, Y internally consistent but not VA)

    # -------------

    # adjust Z for negative VA values in order to have internal consistency, i.e. Y_exports == Y_imports+VA 
    adjustment = ifelse.(VA .< 0.0, VA, 0.0) # spread only negative values
    adjustment = Y .- adjustment
    adjustment = ifelse.(adjustment .< 0.0, 0.0, adjustment)
    adjustment = Y ./ adjustment

    Z = Z .* repeat(adjustment', N*S) # NS×NS, inventory adjustment
    Z = ifelse.(isnan.(Z), 0.0, Z)
    Z = ifelse.(isinf.(Z), 0.0, Z)
    Z = ifelse.(Z .< 0.0, 0.0, Z) # NS×NS

    # rebuild everything
    Y = [sum(Z[i,:]) + sum(F[i,:]) for i in 1:N*S] # NS×1
    VA = [Y[i] - sum(Z[:,i]) for i in 1:N*S] # NS×1, should now be without any negative values

    # calculating gross output via exports/imports+VA should give the same results
    # still not exactly the same! but all differences <1e10
    Y_imports = [sum(Z[:,j]) + VA[j] for j in 1:N*S]
    
    if Y ≈ Y_imports
        println(" ✓ The matrices are internally consistent!")
    else
        println(" × Failure! The matrices are internally inconsistent!")
    end

    
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

    Z_agg = [sum(Z[s:S:(N-1)*S+s,j]) for s in 1:S, j in 1:N*S] # S×NS, sum over origin industries
    π_Z = Z ./ repeat(Z_agg, N) # NS×NS
    π_Z = ifelse.(isnan.(π_Z), 0.0, π_Z) # needed in case Z_agg = 0 (highly unlikely  since all entries in Z >= 0)
    !any(0.0 .<= π_Z .<= 1.0) && println("Problem: not all elements in γ are in intervall [0, 1]")
    
    F_agg = [sum(F[s:S:(N-1)*S+s,j]) for s in 1:S, j in 1:N] # S×N, sum over origin industries
    π_F = F ./ repeat(F_agg, N) # NS×N
    π_F = ifelse.(isnan.(π_F), 0.0, π_F) # needed in case F_agg = 0 (highly unlikely since all entries in F >= 0)
    !any(0.0 .<= π_Z .<= 1.0) && println("Problem: not all elements in γ are in intervall [0, 1]")

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

    # Country-industry level intermediate import expenditure shares
    Z_agg = [sum(Z[i:S:(N-1)*S+i, j]) for i in 1:S, j in 1:N*S] # S×NS, sum over origin industries
    Z_agg_VA = [Z_agg; transpose(VA)] # S+1×NS
    γ_VA = [Z_agg_VA[i,j]/sum(Z_agg_VA[:,j]) for i in 1:S+1, j in 1:N*S]

    γ = γ_VA[1:S, 1:N*S] # S×NS
    γ = ifelse.(isnan.(γ), 0.0, γ) # needed in case Z_agg_VA = 0 (unlikely since all entries in Z, VA >= 0, but still necessary)
    !any(0.0 .<= γ .<= 1.0) && println("Problem: not all elements in γ are in intervall [0, 1]")

    # Country-industry level value added coefficients
    VA_coeff = γ_VA[S+1, 1:N*S] # NS×1
    VA_coeff = ifelse.(isnan.(VA_coeff), 1.0, VA_coeff) # needed in case Z_agg_VA = 0 (unlikely since all entries in Z, VA >= 0, but still necessary)
    #VA_coeff = ifelse.(VA_coeff .< 0.0, 0.0, VA_coeff)

    # Country level value added
    VA_ctry = [sum(VA[i:i+S-1]) for i in 1:S:N*S] # N×1

    # Final imports at country level (sum over rows not columns - otherwise exports)
    F_ctry = [sum(F[:,j]) for j in 1:N] # N×1

    # Country-industry final import expenditure shares (columns sum to 1)
    α = [sum(F[i:S:(N-1)*S+i, j])/F_ctry[j] for i in 1:S, j in 1:N] # S×N
    α = ifelse.(isnan.(α), 0.0, α) # needed in case F_ctry = 0 (unlikely since all entries in F >= 0, but still necessary)
    !any(0.0 .<= α .<= 1.0) && println("Problem: not all elements in γ are in intervall [0, 1]")

    # Country level trade balance (exports - imports)
    # need to calculate country level exports/imports (i.e. aggregate and remove intra-country trade)
    E = [sum(Y[i:i+S-1]) - sum(Z[i:i+S-1,i:i+S-1]) - sum(F[i:i+S-1,ceil(Int, i/S)]) for i in 1:S:N*S] # N×1

    M_Z = [sum(Z[:,j:j+S-1]) - sum(Z[j:j+S-1,j:j+S-1]) for j in 1:S:N*S]
    M_F = [sum(F[:,ceil(Int, j/S)]) - sum(F[j:j+S-1,ceil(Int, j/S)]) for j in 1:S:N*S]
    M = M_Z .+ M_F # N×1

    TB_ctry = E .- M # N×1

    #TB_ctry == VA_ctry .- F_ctry # must hold

    return γ, α, VA_coeff, TB_ctry, VA_ctry, F_ctry
end


# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Summary function
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
    transform_data(dir::String, revision::String, year::Integer, N::Integer, S::Integer)

The function imports and transforms the raw data for further use in the model.

# Arguments
- `dir::String`: directory of raw data (either folder location or entire directory).
- `revision::String`: specifies the WIOD revision (2013 or 2016) as a string.
- `year::Integer`: specifies the year of the WIOD table which should be imported.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

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
julia> Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_data(dir, "2013", 1995, N, S)
```

"""

function transform_data(dir::String, revision::String, year::Integer, N::Integer, S::Integer)

    Z, F, IV = import_data(dir, revision, year, N, S)
    Z, F, Y, VA = inventory_adjustment(Z, F, IV, N, S)
    π_Z, π_F = create_trade_shares(Z, F, N, S)
    γ, α, VA_coeff, TB_ctry, VA_ctry, F_ctry = create_expenditure_shares(Z, F, Y, VA, π_Z, π_F, N, S)

    println(" ✓ Corresponding transformations from WIOD (rev. $revision) for $year were successfully computed! \n")

    return Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F
end

