# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script with functions to compute statistics for estimating sectoral trade elasticities (method of Caliendo and Parro, 2015)
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
    CP_statistic(X::Matrix, τ_X::Matrix, N::Integer, S::Integer)

The function computes the LHS and RHS statistics of equation (23) in Caliendo and Parro (2015).

# Arguments
- `X::Matrix`: NS×N, origin country-industry destination country intermediate/final demand matrix X.
- `τ_X::Matrix`: NS×N, origin country-industry destination country intermediate/final demand tariff matrix τ_X.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `lhs_temp::Vectpr{Float64}`: (Σᴺ⁻²n(n+1)/2)×1, LHS of eq. (23) in logs for intermediate/final demand per industry.
- `rhs_temp::Vectpr{Float64}`: (Σᴺ⁻²n(n+1)/2)×1, RHS of eq. (23) in logs for intermediate/final demand per industry.

# Examples
```julia-repl
julia> lhs_Z_temp, rhs_Z_temp, lhs_F_temp, rhs_F_temp = create_reg_data(X_temp, F_temp, τ_X_temp, τ_F_temp)
```
"""

function CP_statistic(X::Matrix, τ_X::Matrix, N::Integer, S::Integer)

    # initialize loop
    LHS = Float64[]
    RHS = Float64[]

    for i in 1:N-2 # i -> j
        for j in i+1:N-1 # j -> n
            for n in j+1:N # n -> i
                push!( LHS, ( X[i,j]*X[j,n]*X[n,i] ) / ( X[j,i]*X[n,j]*X[i,n] ) )
                push!( RHS, ( τ_X[i,j]*τ_X[j,n]*τ_X[n,i] ) / ( τ_X[j,i]*τ_X[n,j]*τ_X[i,n] ) )
            end
        end
    end

    # create series in logs
    lhs_temp = log.(LHS)
    rhs_temp = log.(RHS)

    return lhs_temp, rhs_temp
end


# ------------


"""
    elasticity_data(M::Matrix, τ::Matrix, demand::String, N::Integer, S::Integer)

A summary function to collect all industry LHS and RHS statistics of equation (23) in Caliendo and Parro (2015). 
    Notice that intermediate and final demand are treated seperately.

# Arguments
- `M::Matrix`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `τ::Matrix`: NS×NS, origin country-industry destination country-industry intermediate input tariff matrix τ_Z.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `lhs::Vectpr{Float64}`: S*(Σᴺ⁻²n(n+1)/2)×1, LHS of eq. (23) in logs for intermediate inputs and all industries.
- `rhs::Vectpr{Float64}`: S*(Σᴺ⁻²n(n+1)/2)×1, RHS of eq. (23) in logs for intermediate inputs and all industries.

# Examples
```julia-repl
julia> lhs, rhs = elasticity_data(Z, F, τ_Z, τ_F, N, S)
```
"""

function elasticity_data(M::Matrix, τ::Matrix, demand::String, N::Integer, S::Integer)

    if demand == "intermediate"
        # reduce intermediate demand to origin country-industry destination country
        X = [sum(M[i,j:j+S-1]) for i in 1:N*S, j in 1:S:N*S] # NS×N

    elseif demand == "total" # wrong
        X = [sum(M[i,j:j+S-1]) for i in 1:N*S, j in 1:S:N*S] # NS×N
        X = X .+ F
    elseif demand == "final" # wrong, dont understand how to reshape final demand matrix
        X = M
    else
        println(" × Failure! Specify demand schedule properly! \n")
    end

    τ_X = τ

    # initialize loop
    lhs = Float64[]
    rhs = Float64[]

    for j in 1:S

        # same sector in rows
        X_temp = X[j:S:(N-1)*S+j, :] # N×N
        τ_X_temp = τ_X[j:S:(N-1)*S+j, :] # N×N
    
        lhs_temp, rhs_temp = CP_statistic(X_temp, τ_X_temp, N, S)

        lhs = [lhs; lhs_temp]
        rhs = [rhs; rhs_temp]
        
    end
    
    # zero trade causes inf => create NaN instead which we let regression filter out later
    lhs = ifelse.(isinf.(lhs), NaN, lhs)
    rhs = ifelse.(isinf.(rhs), NaN, rhs)

    # check if computed length equals number of (theoretical) combinations
    combinations = Int64(sum([n*(n+1)/2 for n in 1:N-2])) # number of possible combinations from formula
    check = ([length(lhs) length(rhs)] .== S*combinations) == [1 1]
    if check == true 
        println(" ✓ The length of the vector coincides with the theoretical length of: $combinations \n")
    else
        println(" × Failure! The length of the vector does not! coincide with the theoretical length of: $combinations \n")
    end

    return lhs, rhs
 
end