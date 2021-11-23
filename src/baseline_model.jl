# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script for the baseline model
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData, LinearAlgebra, Statistics

include("transform_WIOD_2016_3.jl") # Script with functions to import and transform raw data
include("price_hat.jl") # Script with function to obtain the price index
include("wage_hat.jl") # Script with function to obtain the wages and gross output

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

years = 2000:2014 # vector with years covered by WIOD Rev. 2016
N = 44 # number of countries 
S = 56 # number of industries
dir = "C:/Users/u0148308/Desktop/raw/" # location of raw data

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

function baseline(dir::String, year::Integer, N::Integer, S::Integer)
    
    # iteration parameters
    vfactor = 0.4
    tolerance = 1e-3
    max_iteration = 50
    iteration = 0
    max_error = 1e7

    # sectoral trade elasticity
    θ = 5 # assumption
    θ = fill(θ, N*S) # NS×1, work with country-industry elasticities

    # trade costs
    τ_hat_Z = ones(N*S, N*S) # NS×NS
    τ_hat_F = ones(N*S, N) # NS×N

    # initialize wages and price indices
    w_hat = ones(N) # N×1
    #w_hat = fill(2.0, N)
    # ------------

    Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_WIOD_2016(dir, year, N, S)

    #TB_ctry .= 0.0

    # ------------

    while max_error > tolerance && iteration <= max_iteration

        P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat = create_price_index_hat(w_hat, τ_hat_Z, τ_hat_F, VA_coeff, γ, π_Z, π_F, θ)

        # compute counterfactual trade shares
        global π_prime_Z = π_Z .* π_hat_Z # NS×NS
        global π_prime_F = π_F .* π_hat_F # NS×N

        # store last wage in case new optimization obtains negative wages
        w_hat_prev = copy(w_hat)

        w_hat, Y_prime = create_wages_hat(w_hat, vfactor, π_prime_Z, π_prime_F, VA_ctry, TB_ctry, γ, α)

        # update iteration parameters
        error = abs.(w_hat .- w_hat_prev)
        max_error = maximum(error) # update error
        iteration += 1 # update iteration count

        if minimum(w_hat) < 0.0
            w_hat = copy(w_hat_prev)
            println("Iteration $iteration completed with error $max_error (wage negative, rerun with previous estimate)")
        else
            println("Iteration $iteration completed with error $max_error")
        end

    end

    # ------------

    return w_hat, Y_prime, π_prime_Z, π_prime_F, α, γ

end


w_hat, Y_prime, π_prime_Z, π_prime_F, α, γ = baseline(dir, 2014, N, S)
