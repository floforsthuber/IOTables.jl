# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to obtain the wages and gross output
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
    create_wages_hat(w_hat::Vector, vfactor::Number, π_prime_Z::Matrix, π_prime_F::Matrix, VA_ctry::Vector, TB_ctry::Vector, γ::Matrix, α::Matrix)

The function computes the optimal country wage (change) conditional on the economys cost structure by minimizing the distance between successive iterations.

# Arguments
- `w_hat::Vector`: N×1, country wage (change) vector w_hat.
- `vfactor::Number`: scalar ∈ (0,1), scales wage adjustment, the higher the higher the adjustment.
- `π_prime_Z::Matrix`: NS×NS, counterfactual origin country-industry destination country-industry intermediate goods trade share matrix π_Z.
- `π_prime_F::Matrix`: NS×N, counterfactual origin country-industry destination country final goods trade share matrix π_F.
- `VA_ctry::Vector`: N×1, country value added vector VA_ctry.
- `TB_ctry::Vector`: NS×1, country trade balance (X-M) vector TB.
- `γ::Matrix`: S×NS, country-industry intermediate input expenditure share matrix γ.
- `α::Matrix`: S×N, country-industry final expenditure expenditure share matrix α.

# Output
- `P_hat_Z::Matrix{Float64}`: S×NS, country-industry composite intermediate goods price index (change) per industry matrix P_hat_Z.
- `P_hat_F::Matrix{Float64}`: S×N, country-industry composite final goods price index (change) matrix P_hat_F.
- `π_hat_Z::Matrix{Float64}`: NS×NS, origin country-industry destination country-industry intermediate goods trade share (change) matrix π_hat_Z.
- `π_hat_F::Matrix{Float64}`: NS×N, origin country-industry destination country final goods trade share (change) matrix π_hat_F.
- `cost_hat::Matrix{Float64}`: NS×1, country-industry production cost (change) vector cost_hat.

# Examples
```julia-repl
julia> w_hat, Y_prime, ETB_ctry = create_wages_hat(w_hat, vfactor, π_prime_Z, π_prime_F, VA_ctry, TB_ctry, γ, α)
```
"""

function create_wages_hat(w_hat::Vector, vfactor::Number, π_prime_Z::Matrix, π_prime_F::Matrix, VA_ctry::Vector, TB_ctry::Vector, γ::Matrix, α::Matrix)

    F_prime_ctry = w_hat .* VA_ctry .- TB_ctry # N×1, counterfactual country final goods consumption

    # Goods market clearing from equation (35)
    total_sales_F = π_prime_F .* repeat(α,N) .* repeat(transpose(F_prime_ctry),N*S) # NS×N
    total_sales_F = [sum(total_sales_F[i,:]) for i in 1:N*S] # N×1

    A_prime = π_prime_Z .* repeat(γ, N) # NS×NS, intermediate input coefficient matrix
    
    global Y_prime = inv(I - A_prime) * total_sales_F #  NS×1
    Y_prime = ifelse.(Y_prime .< 0.0, 0.0, Y_prime) # Antras and Chor (2018)

    # excess trade balance
    Z_prime = π_prime_Z .* repeat(γ, N) .* repeat(transpose(Y_prime), N*S) # NS×NS
    F_prime = π_prime_F .* repeat(α, N) .* repeat(transpose(F_prime_ctry), N*S) # NS×N

    E_prime = [sum(Y_prime[i:i+S-1]) - sum(Z_prime[i:i+S-1,i:i+S-1]) - sum(F_prime[i:i+S-1,ceil(Int, i/S)]) for i in 1:S:N*S] # N×1

    M_prime_Z = [sum(Z_prime[:,j:j+S-1]) - sum(Z_prime[j:j+S-1,j:j+S-1]) for j in 1:S:N*S]
    M_prime_F = [sum(F_prime[:,ceil(Int, j/S)]) - sum(F_prime[j:j+S-1,ceil(Int, j/S)]) for j in 1:S:N*S]
    M_prime = M_prime_Z .+ M_prime_F # N×1

    global ETB_ctry = E_prime .- M_prime .- TB_ctry # N×1

    # adjust wages to excess trade balance
    norm_ETB_ctry =  ETB_ctry ./ VA_ctry # N×1, normalized excess trade balance
    δ = sign.(norm_ETB_ctry) .* abs.(vfactor.*norm_ETB_ctry) # i.e. if norm_ETB_ctry > 0 => too much exports in model => wages should increase (sign fⁿ gives ± of surplus)
    w_hat = w_hat .* (1.0 .+ δ ./ w_hat) # N×1, increase/decrease wages for countries with an excess surplus/deficit
    
    return w_hat, Y_prime, ETB_ctry
end
