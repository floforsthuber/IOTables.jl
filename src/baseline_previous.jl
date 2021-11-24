# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script to obtain the price index
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData, LinearAlgebra, Statistics

include("transform_WIOD_2016_3.jl") # Script with functions to import and transform raw data

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

years = 2000:2014 # vector with years covered by WIOD Rev. 2016
N = 44 # number of countries 
S = 56 # number of industries
dir = "C:/Users/u0148308/Desktop/raw/" # location of raw data

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_WIOD_2016(dir, 2014, N, S)

# further definitions needed (use same names as Antras and Chor (2018) for the time being)

# iteration parameters
vfactor = 0.4
tolerance = 1e-6
max_iteration = 500
wfmax = 1e7

# sectoral trade elasticity
θ = 5 # assumption
θ = fill(θ, N*S) # NS×1, work with country-industry elasticities

# trade costs
τ_hat_Z = ones(N*S, N*S) # NS×NS
τ_hat_F = ones(N*S, N) # NS×N

# final country level trade balance
# should take this to be either 0 or the final value in 2011?
TB_new = TB_ctry # N×1
TB_ctry .= 0.0

# initialize wages and price indices
w_hat = ones(N) # N×1
#P0_hat = (ones(S,N)./α).^α # S×N
#P0_hat = transpose(P0_hat) # N×S

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Prices
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


function create_price_index_hat(w_hat::Vector, τ_hat_Z::Matrix, τ_hat_F::Matrix, VA_coeff::Vector, γ::Matrix, π_Z::Matrix, π_F::Matrix, θ::Vector)

    # initialize
    w0_hat = w_hat # N×1, update wage vector
    P0_hat_Z = ones(S, N*S) # S×NS, intermediate goods price indices

    # iteration parameters
    tolerance_p = 1e-3
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
        global cost_hat_Z = zeros(N*S, N*S)
        cost_hat = reshape(cost_hat, S, N)
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
        
        cost_Z = π_Z .* cost_hat_Z # NS×NS, origin country-industry destination country-industry price index (inside of summation)

        # sum over origin countries => price index of country-industry composite good in industry
        global P_hat_Z = [sum(cost_Z[i:S:(N-1)*S+i, j])^(-1/θ[i,ceil(Int,j/S)]) for i in 1:S, j in 1:N*S] # S×NS
        P_hat_Z = ifelse.(isinf.(P_hat_Z), 0.0, P_hat_Z) # remove Inf

        # update iteration parameters
        error = abs.(P_hat_Z .- P0_hat_Z)
        max_error_p = maximum(error) # update error
        P0_hat_Z = P_hat_Z # update to new price index
        iteration_p += 1 # update iteration count

        println("Opt. P: Iteration $iteration_p completed with error $max_error_p") # print update on progress

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

    cost_F = π_F .* cost_hat_F # NS×N, the inside of the summation of the price index
    
    # sum over origin countries => price index of country-industry final good
    P_hat_F = [sum(cost_F[i:S:(N-1)*S+i, j])^(-1/θ[i,j]) for i in 1:S, j in 1:N] # S×N
    P_hat_F = ifelse.(isinf.(P_hat_F), 0.0, P_hat_F) # remove Inf

    # ----------
    # trade shares from equation (45) and (46)
    θ = vec(θ)
    θ = repeat(θ, 1, N*S)
    π_hat_Z = (cost_hat_Z ./ repeat(P_hat_Z, N)) .^ (.-θ) # NS×NS

    θ = θ[:, 1:N] # reduce to one NS×N again
    π_hat_F = (cost_hat_F ./ repeat(P_hat_F, N)) .^ (.-θ) # NS×N

    return P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat, cost_hat_Z
end


P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat, cost_hat_Z = create_price_index_hat(w_hat, τ_hat_Z, τ_hat_F, VA_coeff, γ, π_Z, π_F, θ)

# compute counterfactual trade shares
#π_hat_Z = ifelse.(isinf.(π_hat_Z), 0.0, π_hat_Z) # remove Inf
#π_hat_F = ifelse.(isinf.(π_hat_F), 0.0, π_hat_F)

π_prime_Z = π_Z .* π_hat_Z # NS×NS
π_prime_F = π_F .* π_hat_F # NS×N

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Wages, gross output and value added
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

w_hat_prev = copy(w_hat) # store last wage in case new optimization obtains negative wages


function create_wages_hat(w_hat::Vector, vfactor::Number, π_prime_Z::Matrix, π_prime_F::Matrix, VA_ctry::Vector, TB_ctry::Vector, γ::Matrix, α::Matrix)

    F_prime_ctry = w_hat .* VA_ctry .- TB_ctry # N×1, counterfactual country final goods consumption

    # Goods market clearing from equation (35)
    total_sales_F = π_prime_F .* repeat(α,N) .* repeat(transpose(F_prime_ctry),N*S) # NS×N
    total_sales_F = [sum(total_sales_F[i,:]) for i in 1:N*S] # N×1

    A_prime = π_prime_Z .* repeat(γ, N) # NS×NS, intermediate input coefficient matrix
    
    global Y_prime = inv(I - A_prime) * total_sales_F #  NS×1
    Y_prime = ifelse.(Y_prime .< 0.0, 0.0, Y_prime) # Antras and Chor (2018)

    # # excess trade balance
    # Z_prime = π_prime_Z .* repeat(γ, N) .* repeat(transpose(Y_prime), N*S) # NS×NS
    # F_prime = π_prime_F .* repeat(α, N) .* repeat(transpose(F_prime_ctry), N*S) # NS×N

    # E_prime = [sum(Y_prime[i:i+S-1]) - sum(Z_prime[i:i+S-1,i:i+S-1]) - sum(F_prime[i:i+S-1,ceil(Int, i/S)]) for i in 1:S:N*S] # N×1

    # M_prime_Z = [sum(Z_prime[:,j:j+S-1]) - sum(Z_prime[j:j+S-1,j:j+S-1]) for j in 1:S:N*S]
    # M_prime_F = [sum(F_prime[:,ceil(Int, j/S)]) - sum(F_prime[j:j+S-1,ceil(Int, j/S)]) for j in 1:S:N*S]
    # M_prime = M_prime_Z .+ M_prime_F # N×1

    # ETB_ctry = E_prime .- M_prime .- TB_ctry # N×1

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

    ETB_ctry = RHS .- LHS .- TB_ctry # N×1, excess trade balance

    # # adjust wages to excess trade balance
    # norm_ETB_ctry =  ETB_ctry ./ VA_ctry # N×1, normalized excess trade balance
    # δ = sign.(norm_ETB_ctry) .* abs.(vfactor.*norm_ETB_ctry) # i.e. if norm_ETB_ctry > 0 => too much exports in model => wages should increase (sign fⁿ gives ± of surplus)
    # w_hat = w_hat .* (1.0 .+ δ ./ w_hat) # N×1, increase/decrease wages for countries with an excess surplus/deficit


    # Problem: After first wage adjustment rounds the counterfactual output explodes (Y_prime)
    # which causes the trade balance to be huge and VA_ctry is not enough to force: norm_ETB_ctry ∈ [0,1]
    # hence wages increase by crazy amount
    # Solution: normalize by using Y_prime_ctry instead of VA_ctry, could also recalculate VA_prime_ctry?

    # Problem: at some point the scaling of the adjustment δ ./ w_hat > 1 if w_hat is small enough
    # hence w_hat becomes negative and we run in an infinitive loop
    # Solution: do not scale adjustment by w_hat

    # adjust wages to excess trade balance
    Y_prime_ctry = [sum(Y_prime[i:i+S-1]) for i in 1:S:N*S] # N×1
    norm_ETB_ctry =  ETB_ctry ./ Y_prime_ctry # N×1, normalized excess trade balance
    δ = sign.(norm_ETB_ctry) .* abs.(vfactor.*norm_ETB_ctry) # i.e. if norm_ETB_ctry > 0 => too much exports in model => wages should increase (sign fⁿ gives ± of surplus)
    w_hat = w_hat .* (1.0 .+ δ) # N×1, increase/decrease wages for countries with an excess surplus/deficit

    return w_hat, Y_prime
end

w_hat, Y_prime = create_wages_hat(w_hat, vfactor, π_prime_Z, π_prime_F, VA_ctry, TB_ctry, γ, α)


# second time 
P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat = create_price_index_hat(w_hat, τ_hat_Z, τ_hat_F, VA_coeff, γ, π_Z, π_F, θ)

π_prime_Z = π_Z .* π_hat_Z # NS×NS
π_prime_F = π_F .* π_hat_F # NS×N
w_hat_prev = copy(w_hat)

w_hat, Y_prime = create_wages_hat(w_hat, vfactor, π_prime_Z, π_prime_F, VA_ctry, TB_ctry, γ, α)

# third time 

P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat = create_price_index_hat(w_hat, τ_hat_Z, τ_hat_F, VA_coeff, γ, π_Z, π_F, θ)

π_prime_Z = π_Z .* π_hat_Z # NS×NS
π_prime_F = π_F .* π_hat_F # NS×N
w_hat_prev = copy(w_hat)

w_hat, Y_prime = create_wages_hat(w_hat, vfactor, π_prime_Z, π_prime_F, VA_ctry, TB_ctry, γ, α)
print(findall(w_hat.<0))

# fourth time 

P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat = create_price_index_hat(w_hat, τ_hat_Z, τ_hat_F, VA_coeff, γ, π_Z, π_F, θ)

π_prime_Z = π_Z .* π_hat_Z # NS×NS
π_prime_F = π_F .* π_hat_F # NS×N
w_hat_prev = copy(w_hat)

w_hat, Y_prime = create_wages_hat(w_hat, vfactor, π_prime_Z, π_prime_F, VA_ctry, TB_ctry, γ, α)
print(findall(w_hat.<0))

# fith time 

P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat = create_price_index_hat(w_hat, τ_hat_Z, τ_hat_F, VA_coeff, γ, π_Z, π_F, θ)

π_prime_Z = π_Z .* π_hat_Z # NS×NS
π_prime_F = π_F .* π_hat_F # NS×N
w_hat_prev = copy(w_hat)

w_hat, Y_prime = create_wages_hat(w_hat, vfactor, π_prime_Z, π_prime_F, VA_ctry, TB_ctry, γ, α)
print(findall(w_hat.<0))

### -------------

F_prime_ctry = w_hat_prev .* VA_ctry .- TB_ctry # N×1, counterfactual country final goods consumption

# Goods market clearing from equation (35)
total_sales_F = π_prime_F .* repeat(α,N) .* repeat(transpose(F_prime_ctry),N*S) # NS×N
total_sales_F = [sum(total_sales_F[i,:]) for i in 1:N*S] # N×1

A_prime = π_prime_Z .* repeat(γ, N) # NS×NS, intermediate input coefficient matrix

Y_prime = inv(I - A_prime) * total_sales_F #  NS×1
Y_prime = ifelse.(Y_prime .< 0.0, 0.0, Y_prime) # Antras and Chor (2018)

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

ETB_ctry = RHS .- LHS .- TB_ctry # N×1, excess trade balance

Y_prime_ctry = [sum(Y_prime[i:i+S-1]) for i in 1:S:N*S] # N×1
norm_ETB_ctry =  ETB_ctry ./ (Y_prime_ctry .- ETB_ctry) # N×1, normalized excess trade balance
δ = sign.(norm_ETB_ctry) .* abs.(vfactor.*norm_ETB_ctry) # i.e. if norm_ETB_ctry > 0 => too much exports in model => wages should increase (sign fⁿ gives ± of surplus)
w_hat = w_hat .* (1.0 .+ δ) # N×1, increase/decrease wages for countries with an excess surplus/deficit

ETB_ctry[9]
Y_prime_ctry[9]

sum(Y_prime[S*9:S*10-1])