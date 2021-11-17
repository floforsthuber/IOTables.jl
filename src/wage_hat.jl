# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to obtain the wages and gross output
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


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

    #### Antras and Chor (2018) calculation
    LHS = zeros(N)
    RHS = zeros(N)
    for j in 1:N
        for s in 1:N
            for i in 1:N
                iiss = (i-1)*S+s
                jjss = (j-1)*S+s
                for r in 1:S
                    jjrr = (j-1)*S+r
                    iirr = (i-1)*S+r
                    LHS[j] += π_prime_Z[iiss,jjrr]*γ[s,jjrr]*Y_prime[jjrr]
                    RHS[j] += π_prime_Z[jjss,iirr]*γ[s,iirr]*Y_prime[iirr]
                end
                LHS[j] += π_prime_F[iiss,j]*α[s,j]*F_prime_ctry[j]
                RHS[j] += π_prime_F[jjss,i]*α[s,i]*F_prime_ctry[i]
            end
        end
    end

    ETB_ctry = E_prime .- M_prime .- TB_ctry # N×1
    ETB_ctry = RHS .- LHS .- TB_ctry # N×1, excess trade balance

    # adjust wages to excess trade balance
    norm_ETB_ctry =  ETB_ctry ./ VA_ctry # N×1, normalized excess trade balance
    δ = sign.(norm_ETB_ctry) ./ abs.(vfactor.*norm_ETB_ctry) # i.e. if norm_ETB_ctry > 0 => too much exports in model => wages should increase (sign fⁿ gives ± of surplus)
    w_hat = w_hat .* (1.0 .+ δ ./ w_hat) # N×1, increase/decrease wages for countries with an excess surplus/deficit

    return w_hat, Y_prime
end
