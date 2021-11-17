# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to obtain the price index
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
    create_price_index_hat(w_hat::Vector, τ_hat_Z::Matrix, τ_hat_F::Matrix, VA_coeff::Vector, γ::Matrix, π_Z::Matrix, π_F::Matrix, θ::Vector)

The function performs computes the optimal price indices conditional on the economy's wage structure by minimizing the distance between successive iterations.

# Arguments
- `w_hat::Vector`: N×1, country wages vector w_hat.
- `τ_hat_Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate goods demand trade cost matrix τ_hat_Z.
- `τ_hat_F::Matrix`: NS×N, origin country-industry destination country final goods trade cost matrix τ_hat_F.
- `VA_coeff::Vector`: NS×1, country-industry value added coefficients vector VA_coeff.
- `γ::Matrix`: S×NS, country-industry intermediate input expenditure share matrix γ.
- `π_Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate goods trade share matrix π_Z.
- `π_F::Matrix`: NS×N, origin country-industry destination country final goods trade share matrix π_F.
- `θ::Vector`: NS×1, country-industry trade elasticities vector θ.

# Output
- `P_hat_Z::Matrix{Float64}`: S×NS, country-industry composite intermediate goods price index (change) per industry matrix P_hat_Z.
- `P_hat_F::Matrix{Float64}`: S×N, country-industry composite final goods price index (change) matrix P_hat_F.
- `π_hat_Z::Matrix{Float64}`: NS×NS, origin country-industry destination country-industry intermediate goods trade share (change) matrix π_hat_Z.
- `π_hat_F::Matrix{Float64}`: NS×N, origin country-industry destination country final goods trade share (change) matrix π_hat_F.
- `cost_hat::Matrix{Float64}`: NS×1, country-industry production cost (change) vector cost_hat.

# Examples
```julia-repl
julia> P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat = create_price_index_hat(w_hat, τ_hat_Z, τ_hat_F, VA_coeff, γ, π_Z, π_F, θ)
```
"""

function create_price_index_hat(w_hat::Vector, τ_hat_Z::Matrix, τ_hat_F::Matrix, VA_coeff::Vector, γ::Matrix, π_Z::Matrix, π_F::Matrix, θ::Vector)

    # initialize
    w0_hat = w_hat # N×1, update wage vector
    P0_hat_Z = ones(S, N*S) # S×NS, intermediate goods price indices

    # iteration parameters
    tolerance_p = 1e-1
    max_iteration_p = 50
    iteration_p = 0
    max_error_p = 1.0

    while max_error_p > tolerance_p && iteration_p <= max_iteration_p

        # cost structure of the economy, equation (47)
        cost_hat_w = [w0_hat[ceil(Int,i/S)]^VA_coeff[i] for i in 1:N*S] # NS×1, wages
        cost_hat_p = P0_hat_Z.^γ
        cost_hat_p = [prod(cost_hat_p[:,i]) for i in 1:N*S] # NS×1, price indices (prod is the same as sum just for multiplication)
        global cost_hat = cost_hat_w .* cost_hat_p # NS×1

        # ----------
        # price index for intermediate goods from equation (48)
        global cost_Z = [π_Z[i,j]*(cost_hat[i]*τ_hat_Z[i,j])^(-θ[i]) for i in 1:N*S, j in 1:N*S] # NS×NS

        # sum over origin countries => price index of country-industry composite good in industry
        global cost_agg_Z = zeros(S, N*S) # S×NS
        global P_hat_Z = zeros(S, N*S) # S×NS
        θ = reshape(θ, S, N) # S×N
        for i in 1:S
            for j in 1:N*S
                for k in 0:S:N*S-1
                    cost_agg_Z[i, j] += cost_Z[k+i,j]
                end
                P_hat_Z[i, j] = cost_agg_Z[i, j]^(-1/θ[i,ceil(Int,j/S)])
            end
        end

        # update iteration parameters
        error = abs.(P_hat_Z .- P0_hat_Z)
        max_error_p = maximum(error) # update error
        P0_hat_Z = P_hat_Z # update to new price index
        iteration_p += 1 # update iteration count

        println("Opt. P: Iteration $iteration_p completed with error $max_error_p") # print update on progress

    end


    # ----------
    # price index for final goods from equation (49) --- not used in optimization, i.e. residual optimum (since cost is same in F as in Z!)
    θ = vec(θ) # NS×1, reshape to have row vector
    cost_F = [π_F[i,j]*(cost_hat[i]*τ_hat_F[i,j])^(-θ[i]) for i in 1:N*S, j in 1:N] # NS×N

    cost_agg_F = zeros(S, N) # S×N
    P_hat_F = zeros(S, N) # S×N
    θ = reshape(θ, S, N) # S×N
    for i in 1:S
        for j in 1:N
            for k in 0:S:N*S-1
                cost_agg_F[i, j] += cost_F[k+i,j]
            end
            P_hat_F[i,j] = cost_agg_F[i, j]^(-1/θ[i,j])
        end
    end


    # ----------
    # trade shares from equation (45) and (46)
    π_hat_Z = cost_Z./repeat(cost_agg_Z, N) # NS×NS
    π_hat_F = cost_F./repeat(cost_agg_F, N) # NS×N

    π_hat_Z = ifelse.(isinf.(π_hat_Z), 0.0, π_hat_Z) # remove Inf
    π_hat_F = ifelse.(isinf.(π_hat_F), 0.0, π_hat_F)


    return P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat
end
