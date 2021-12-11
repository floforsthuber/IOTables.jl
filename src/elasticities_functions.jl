# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script with functions to compute statistics for estimating sectoral trade elasticities (method of Caliendo and Parro, 2015)
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
    CP_statistic(X::Matrix, F::Matrix, τ_X::Matrix, τ_F::Matrix, N::Integer, S::Integer)

The function computes the LHS and RHS statistics of equation (23) in Caliendo and Parro (2015). Notice that intermediate and final demand 
    are treated seperately with the goal in mind to obtain seperate trade elasticities for intermediate and final demand.

# Arguments
- `X::Matrix`: NS×N, origin country-industry destination country intermediate demand matrix Z.
- `F::Matrix`: NS×N, origin country-industry destination country final demand matrix F.
- `τ_Z::Matrix`: NS×N, origin country-industry destination country intermediate input tariff matrix τ_Z.
- `τ_F::Matrix`: NS×NS, origin country-industry destination country final demand tariff matrix τ_F.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `lhs_Z_temp::Vectpr{Float64}`: (Σᴺ⁻²n(n+1)/2)×1, LHS of eq. (23) in logs for intermediate inputs per industry.
- `rhs_Z_temp::Vectpr{Float64}`: (Σᴺ⁻²n(n+1)/2)×1, RHS of eq. (23) in logs for intermediate inputs per industry.
- `lhs_F_temp::Vectpr{Float64}`: (Σᴺ⁻²n(n+1)/2)×1, LHS of eq. (23) in logs for final demand per industry.
- `rhs_F_temp::Vectpr{Float64}`: (Σᴺ⁻²n(n+1)/2)×1, RHS of eq. (23) in logs for final demand per industry.

# Examples
```julia-repl
julia> lhs_Z_temp, rhs_Z_temp, lhs_F_temp, rhs_F_temp = create_reg_data(X_temp, F_temp, τ_X_temp, τ_F_temp)
```
"""

function CP_statistic(X::Matrix, F::Matrix, τ_X::Matrix, τ_F::Matrix, N::Integer, S::Integer)

    # initialize loop
    LHS_Z = Float64[]
    RHS_Z = Float64[]

    LHS_F = Float64[]
    RHS_F = Float64[]

    for i in 1:N-2 # i -> j
        for j in i+1:N-1 # j -> n
            for n in j+1:N # n -> i
                push!( LHS_Z, ( X[i,j]*X[j,n]*X[n,i] ) / ( X[j,i]*X[n,j]*X[i,n] ) )
                push!( RHS_Z, ( τ_X[i,j]*τ_X[j,n]*τ_X[n,i] ) / ( τ_X[j,i]*τ_X[n,j]*τ_X[i,n] ) )

                push!( LHS_F, ( F[i,j]*F[j,n]*F[n,i] ) / ( F[j,i]*F[n,j]*F[i,n] ) )
                push!( RHS_F, ( τ_F[i,j]*τ_F[j,n]*τ_F[n,i] ) / ( τ_F[j,i]*τ_F[n,j]*τ_F[i,n] ) )
            end
        end
    end

    # create series in logs
    lhs_Z_temp = log.(LHS_Z)
    rhs_Z_temp = log.(RHS_Z)

    lhs_F_temp = log.(LHS_F)
    rhs_F_temp = log.(RHS_F)


    return lhs_Z_temp, rhs_Z_temp, lhs_F_temp, rhs_F_temp
end


# ------------


"""
    CP_statistic(Z::Matrix, F::Matrix, τ_Z::Matrix, τ_F::Matrix, N::Integer, S::Integer)

A summary function to collect all industry LHS and RHS statistics of equation (23) in Caliendo and Parro (2015). Notice that intermediate and final demand 
    are treated seperately.

# Arguments
- `Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: NS×N, origin country-industry destination country final demand matrix F.
- `τ_Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate input tariff matrix τ_Z.
- `τ_F::Matrix`: NS×NS, origin country-industry destination country final demand tariff matrix τ_F.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `lhs_Z::Vectpr{Float64}`: S*(Σᴺ⁻²n(n+1)/2)×1, LHS of eq. (23) in logs for intermediate inputs and all industries.
- `rhs_Z::Vectpr{Float64}`: S*(Σᴺ⁻²n(n+1)/2)×1, RHS of eq. (23) in logs for intermediate inputs and all industries.
- `lhs_F::Vectpr{Float64}`: S*(Σᴺ⁻²n(n+1)/2)×1, LHS of eq. (23) in logs for final demand and all industries.
- `rhs_F::Vectpr{Float64}`: S*(Σᴺ⁻²n(n+1)/2)×1, RHS of eq. (23) in logs for final demand and all industries.

# Examples
```julia-repl
julia> lhs_Z, rhs_Z, lhs_F, rhs_F = elasticity_data(Z, F, τ_Z, τ_F, N, S)
```
"""

function elasticity_data(Z::Matrix, F::Matrix, τ_Z::Matrix, τ_F::Matrix, N::Integer, S::Integer)

    # reduce intermediate demand to origin country-industry destination country
    X = [sum(Z[i,j:j+S-1]) for i in 1:N*S, j in 1:S:N*S] # NS×N
    #X = X .+ F

    if size(τ_Z, 2) == N*S
        τ_X = [mean(τ_Z[i,j:j+S-1]) for i in 1:N*S, j in 1:S:N*S] # NS×N
    else
        τ_X = τ_Z
    end

    # initialize loop
    lhs_Z = Float64[]
    rhs_Z = Float64[]

    lhs_F = Float64[]
    rhs_F = Float64[]

    for j in 1:S

        # same sector in rows
        X_temp = X[j:S:(N-1)*S+j, :] # S×N
        τ_X_temp = τ_X[j:S:(N-1)*S+j, :] # S×N
        F_temp = F[j:S:(N-1)*S+j, :] # S×N
        τ_F_temp = τ_F[j:S:(N-1)*S+j, :] # S×N
    
        lhs_Z_temp, rhs_Z_temp, lhs_F_temp, rhs_F_temp = CP_statistic(X_temp, F_temp, τ_X_temp, τ_F_temp, N, S)

        lhs_Z = [lhs_Z; lhs_Z_temp]
        rhs_Z = [rhs_Z; rhs_Z_temp]
        lhs_F = [lhs_F; lhs_F_temp]
        rhs_F = [rhs_F; rhs_F_temp]
        
    end
    
    # zero trade causes inf => create NaN instead which we let regression filter out later
    lhs_Z = ifelse.(isinf.(lhs_Z), NaN, lhs_Z)
    rhs_Z = ifelse.(isinf.(rhs_Z), NaN, rhs_Z)
    
    lhs_F = ifelse.(isinf.(lhs_F), NaN, lhs_F)
    rhs_F = ifelse.(isinf.(rhs_F), NaN, rhs_F)

    # check if computed length equals number of (theoretical) combinations
    combinations = Int64(sum([n*(n+1)/2 for n in 1:N-2])) # number of possible combinations from formula
    check = ([length(lhs_Z) length(rhs_Z) length(lhs_F) length(rhs_F)] .== S*combinations) == [1 1 1 1]
    if check == true 
        println(" ✓ The length of the vector coincides with the theoretical length of: $combinations \n")
    else
        println(" × Failure! The length of the vector does not! coincide with the theoretical length of: $combinations \n")
    end

    return lhs_Z, rhs_Z, lhs_F, rhs_F
 
end