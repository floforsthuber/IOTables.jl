
# --------------- Subsetting IO table -----------------------------------------------------------------------------------------------------------------------------

"""
    index_subset(N::Integer, S::Integer, ctry_all::Vector{String}, ctry_subset::Vector{String}, ind_all::Vector{String}, ind_subset::Vector{String})

The function yields the row/column indices for the IO table according to desired countries and industries. The subsetting should be done for Z and F,
    which then needs a recalculation of Y and VA!

# Arguments
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.
- `ctry_all::Vector{String}`: ordered vector of all country names in IO table.
- `ctry_subset::Vector{String}`: vector of country names to subset.
- `ind_all::Vector{String}`: ordered vector of all industry codes (2 digits) in IO table.
- `ind_subset::Vector{String}`: vector of industry codes (2 digits) to subset.

# Output
- `index_subset::Vector{Integer}`: indices to subset Z and F matrices.

# Examples
```julia-repl
julia> index_subset = subset(N, S, countries, ctry_EU, industries, ind_sub)
```
"""

function index_subset(N::Integer, S::Integer, ctry_all::Vector{String}, ctry_subset::Vector{String}, ind_all::Vector{String}, ind_subset::Vector{String})
    
    # subsetting rows
    N_subset = length(ctry_subset)
    S_subset = length(ind_subset)

    all = repeat(ctry_all, inner=S) .* "_" .* repeat(ind_all, outer=N) # all row country-industry pairs
    subset = repeat(ctry_subset, inner=S_subset) .* "_" .* repeat(ind_subset, outer=N_subset) # subset of country-industry pairs
    index_subset = findall(in(subset), all)

    return index_subset
end


# ---------------

"""
    ctry_aggregation(Z::Matrix, F::Matrix, N::Integer, S::Integer, 
        ctry_all::Vector{String}, ctry_subset::Vector{String}, ind_all::Vector{String}, ind_subset::Vector{String})

The function creates new intermediate and final demand matrices by aggregating the countries in the vector "ctry_subset". The aggregation only works for 
    countries, therefore, the two industry vectors need to be the same.

# Arguments
- `Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: NS×N, origin country-industry destination country final demand matrix F.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.
- `ctry_all::Vector{String}`: ordered vector of all country names in IO table.
- `ctry_subset::Vector{String}`: vector of country names to subset.
- `ind_all::Vector{String}`: ordered vector of all industry codes (2 digits) in IO table.
- `ind_subset::Vector{String}`: needs to be identical to the vector "ind_all".

# Output
- `Z::Matrix`: N_new*S×N_new*S, new origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: N_new*S×N_new, new origin country-industry destination country final demand matrix F.

# Examples
```julia-repl
julia> Z_new, F_new = aggregation(Z, F, N, S, countries, ctry_EU, industries, industries)
```
"""

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





# ---------------- aggregation in one function --------------------------------------------------------------------------------------------------------------------------

"""
    ctry_aggregation(Z::Matrix, F::Matrix, N::Integer, S::Integer, ctry_all::Vector{String}, ctry_subset::Vector{String})

The function creates new intermediate and final demand matrices by aggregating the countries in the vector "ctry_subset". The aggregation only works for 
    countries not for industries.

# Arguments
- `Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: NS×N, origin country-industry destination country final demand matrix F.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.
- `ctry_all::Vector{String}`: ordered vector of all country names in IO table.
- `ctry_subset::Vector{String}`: vector of country names to subset.

# Output
- `Z::Matrix`: N_new*S×N_new*S, new origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: N_new*S×N_new, new origin country-industry destination country final demand matrix F.

# Examples
```julia-repl
julia> Z_new, F_new = aggregation(Z, F, N, S, countries, ctry_EU)
```
"""

function ctry_aggregation(Z::Matrix, F::Matrix, N::Integer, S::Integer, ctry_all::Vector{String}, ctry_subset::Vector{String})

    # obtain indices for subsetting Z and F
    N_subset = length(ctry_subset)
    industries = lpad.(1:S, 2, '0') # vector of industries

    all = repeat(ctry_all, inner=S) .* "_" .* repeat(industries, outer=N) # all row country-industry pairs
    subset = repeat(ctry_subset, inner=S) .* "_" .* repeat(industries, outer=N_subset) # subset of country-industry pairs
    
    index_subset = findall(in(subset), all) # can be used for both rows and columns since IO table is a square matrix
    index_subset_F = findall(in(ctry_subset), ctry_all) # only have countries in columns

    # subset columns, construct aggregate and rebuild matrix
    Z_col_subset = Z[:, index_subset]
    Z_col_subset = [sum(Z_col_subset[i, j:S:(N_subset-1)*S+j]) for i in 1:N*S, j in 1:S] # sum over industries to create aggregate
    Z_col_other = Z[:, Not(index_subset)] # other countries
    Z_new = hcat(Z_col_other, Z_col_subset) # rebuild Z with new aggregate as last columns

    F_col_subset = F[:, index_subset_F]
    F_col_subset = [sum(F_col_subset[i,:]) for i in 1:N*S]
    F_col_other = F[:, Not(index_subset_F)] # other countries
    F_new = hcat(F_col_other, F_col_subset) # rebuild F with new aggregate as last column

    # subset rows from new matrix, construct aggregate and rebuild matrix
    N_new = size(F_new, 2) # new number of countries

    Z_row_subset = Z_new[index_subset, :]
    Z_row_subset = [sum(Z_row_subset[i:S:(N_subset-1)*S+i, j]) for i in 1:S, j in 1:N_new*S]
    Z_row_other = Z_new[Not(index_subset), :]
    Z_new = vcat(Z_row_other, Z_row_subset)

    F_row_subset = F_new[index_subset, :]
    F_row_subset = [sum(F_row_subset[i:S:(N_subset-1)*S+i, j]) for i in 1:S, j in 1:N_new]
    F_row_other = F_new[Not(index_subset), :]
    F_new = vcat(F_row_other, F_row_subset)

    return Z_new, F_new
end