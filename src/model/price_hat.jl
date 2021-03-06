# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to obtain the price index
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
    create_price_index_hat(w_hat::Vector, τ_hat_Z::Matrix, τ_hat_F::Matrix, VA_coeff::Vector, γ::Matrix, π_Z::Matrix, π_F::Matrix, θ::Vector)

The function computes the optimal price indices (change) conditional on the economys' wage structure by minimizing the distance between successive iterations 
    in the intermediate input price index.

# Arguments
- `w_hat::Vector`: N×1, country wage (change) vector w_hat.
- `τ_hat_Z::Matrix`: NS×NS, origin country-industry destination country-industry intermediate goods demand trade cost (change) matrix τ_hat_Z.
- `τ_hat_F::Matrix`: NS×N, origin country-industry destination country final goods trade cost (change) matrix τ_hat_F.
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
    w0_hat = copy(w_hat) # N×1, update wage vector
    P0_hat_Z = ones(S, N*S) # S×NS, intermediate goods price indices

    # iteration parameters
    tolerance_p = 1e-3
    max_iteration_p = 100
    iteration_p = 0
    max_error_p = 1.0

    while max_error_p > tolerance_p && iteration_p <= max_iteration_p

        # # cost structure of the economy, equation (47)
        # cost_hat_w = [w0_hat[ceil(Int,i/S)]^VA_coeff[i] for i in 1:N*S] # NS×1, wages
        # cost_hat_p = P0_hat_Z.^γ
        # cost_hat_p = [prod(cost_hat_p[:,i]) for i in 1:N*S] # NS×1, price indices (prod is the same as sum just for multiplication)
        # cost_hat = cost_hat_w .* cost_hat_p # NS×1
        # global cost_hat = reshape(cost_hat, S, N) # S×N

        global cost_hat = ones(S,N)
        VAcoeff = reshape(VA_coeff, S, N)
        for s in 1:S
            for j in 1:N
                cost_hat[s,j] = w0_hat[j]^VAcoeff[s,j]
                for r in 1:S
                    cost_hat[s,j] = cost_hat[s,j]*P0_hat_Z[r,(j-1)*S+s]^γ[r,(j-1)*S+s]
                end
            end
        end

        # ----------
        # price index for intermediate goods from equation (48)
        global cost_hat_Z = zeros(N*S, N*S)
        θ = reshape(θ, S, N) # S×N

        for i in 1:N
            for r in 1:S
                for j in 1:N
                    for s in 1:S
                        iirr = (i-1)*S+r
                        jjss = (j-1)*S+s
                        cost_hat_Z[iirr,jjss] = (cost_hat[r,i]*τ_hat_Z[iirr,jjss])^(-θ[r,i])
                    end
                end
            end
        end

        # price index calculation
        global P_hat_Z = zeros(S, N*S)
        for r in 1:S
            for j in 1:N
                for s in 1:S
                    jjss = (j-1)*S+s
                    P_hat_Z[r,jjss] = (sum(π_Z[r:S:(N-1)*S+r,jjss] .* cost_hat_Z[r:S:(N-1)*S+r,jjss]))^(-1/θ[r,j])
                end
            end
        end

        P_hat_Z = ifelse.(isinf.(P_hat_Z), 0.0, P_hat_Z)

        # update iteration parameters
        error = abs.(P_hat_Z .- P0_hat_Z)
        max_error_p = maximum(error) # update error
        P0_hat_Z[:] = copy(P_hat_Z) # update to new price index
        iteration_p += 1 # update iteration count

        println(" - Inner loop: Iteration $iteration_p completed with error $max_error_p") # print update on progress

    end

    # ----------
    # price index for final goods from equation (49) --- not used in optimization, i.e. residual optimum (since cost is same in F as in Z!)
    cost_hat_F = zeros(N*S,N)
    for i in 1:N
        for r in 1:S
            for j in 1:N
                iirr = (i-1)*S+r
                cost_hat_F[iirr,j] = (cost_hat[r,i]*τ_hat_F[iirr,j])^(-θ[r,i])
            end
        end
    end

    global P_hat_F = zeros(S, N)
    for r in 1:S
        for j in 1:N
            P_hat_F[r,j] = (sum(π_F[r:S:(N-1)*S+r,j] .* cost_hat_F[r:S:(N-1)*S+r,j]))^(-1/θ[r,j])
        end
    end

    P_hat_F = ifelse.(isinf.(P_hat_F), 0.0, P_hat_F)

    # ----------
    # trade shares from equation (45) and (46)
    θ = vec(θ)
    θ = repeat(θ, 1, N*S)
    global π_hat_Z = cost_hat_Z ./ repeat(P_hat_Z, N) .^ (.-θ) # NS×NS, - cost_hat_Z is already to power of -θ => only P_hat_Z^-θ

    θ = θ[:, 1:N] # reduce to one NS×N again
    global π_hat_F = cost_hat_F ./ repeat(P_hat_F, N) .^ (.-θ) # NS×N,  - cost_hat_Z is already to power of -θ => only P_hat_Z^-θ

    return P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat
end
