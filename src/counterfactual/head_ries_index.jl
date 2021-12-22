# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script with functions to create bilateral Head-Ries index (symmetric bilateral trade costs)
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
    head_ries_index(Z::Matrix, F::Matrix, θ, N::Integer, S::Integer)

The function computes bilateral trade costs according to Head-Ries, i.e. symmetric inter country-industry trade costs and 
    zero intra country trade costs. Trade costs are computed for both intermediate input and final demand. The trade elasticity can 
    take the form of a uniform country-industry value or origin country-industry trade elasticity.

# Arguments
- `Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: NS×N, origin country-industry destination country final demand matrix F.
- `θ`: assumption for the effective trade elasticity, can either be supplied as a simple value or a S×N matrix.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `τ_HR_Z::Matrix{Float64}`: NS×NS, Head-Reis Index, symmetric bilateral country-industry intermediate input trade costs 
        with intra-country trade assumed to have no trade costs.
- `τ_HR_F::Matrix{Float64}`: NS×N,  Head-Reis Index, symmetric bilateral country-industry final demand trade costs 
    with intra-country trade assumed to have no trade costs.


# Examples
```julia-repl
julia> τ_HR_Z, τ_HR_F = head_ries_index(Z, F, θ, N, S)
```
"""

function head_ries_index(Z::Matrix, F::Matrix, θ, N::Integer, S::Integer)

    if length(θ) == 1
        θ = fill(θ, N*S) # NS×1, work with country-industry elasticities
        θ = reshape(θ, S, N) # S×N
    elseif size(θ) == (S, N)
    else
        println("The effective trade elasticity θ is not supplied in the right format: either a simple number or a S×N matrix.")
    end

    Z_HR = ifelse.(Z .== 0.0, 1e-18, Z) # otherwise we would have infinite trade costs
    τ_HR_Z = zeros(N*S, N*S) # NS×NS

    for ctry1 in 1:N
        for ctry2 in 1:N
            for ind1 in 1:S
                for ind2 in 1:S
                    if ctry1 == ctry2
                        τ_HR_Z[ind1+(ctry1-1)*S, ind2+(ctry2-1)*S] = 1.0
                    else
                        Z1122 = Z_HR[ind1+(ctry1-1)*S, ind2+(ctry2-1)*S]
                        Z2122 = Z_HR[ind1+(ctry2-1)*S, ind2+(ctry2-1)*S]
                        Z2112 = Z_HR[ind1+(ctry2-1)*S, ind2+(ctry1-1)*S]
                        Z1112 = Z_HR[ind1+(ctry1-1)*S, ind2+(ctry1-1)*S]
                        τ_temp = ((Z1122*Z2112)/(Z2122*Z1112))^(-1/(2*θ[ind1,ctry1]))
                        τ_HR_Z[ind1+(ctry1-1)*S, ind2+(ctry2-1)*S] = τ_temp
                        τ_HR_Z[ind1+(ctry2-1)*S, ind2+(ctry1-1)*S] = τ_temp
                    end
                end
            end
        end
    end

    F_HR = ifelse.(F .== 0.0, 1e-18, F)
    τ_HR_F = zeros(N*S, N) # NS×N

    for ctry1 in 1:N
        for ctry2 in 1:N
            for ind1 in 1:S
                if ctry1 == ctry2
                    τ_HR_F[ind1+(ctry1-1)*S, ctry2] = 1.0
                else
                    F12 = F_HR[ind1+(ctry1-1)*S, ctry2]
                    F22 = F_HR[ind1+(ctry2-1)*S, ctry2]
                    F21 = F_HR[ind1+(ctry2-1)*S, ctry1]
                    F11 = F_HR[ind1+(ctry1-1)*S, ctry1]
                    τ_temp = ((F12*F21)/(F11*F22))^(-1/(2*θ[ind1,ctry1]))
                    τ_HR_F[ind1+(ctry1-1)*S, ctry2] = τ_temp
                    τ_HR_F[ind1+(ctry2-1)*S, ctry1] = τ_temp
                end
            end
        end
    end

    return τ_HR_Z, τ_HR_F
end
