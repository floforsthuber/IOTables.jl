# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Baseline model
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData,  XLSX, LinearAlgebra, Statistics

dir = "X:/VIVES/1-Personal/Florian/git/IOTables/src/"
include(dir * "model/transform_data.jl") # Script with functions to import and transform raw data
include(dir * "model/price_hat.jl") # Script with function to obtain the price index
include(dir * "model/wage_hat.jl") # Script with function to obtain the wages and gross output

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

dir = "C:/Users/u0148308/data/raw/" # location of raw data

# Data specification

# WIOD rev. 2013
source = "WIOD"
revision = "2013"
year = 1995 # specified year
N = 41 # number of countries 
S = 35 # number of industries

# # WIOD rev. 2016
# source = "WIOD"
# revision = "2016"
# year = 2014 # specified year
# N = 44 # number of countries 
# S = 56 # number of industries


# -------------------------------------------------------------------------------------------------------------------------------------------------------------

Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_data(dir, source, revision, year, N, S)

# iteration parameters
vfactor = 0.4
tolerance = 1e-6
max_iteration = 100
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

# # adjust trade balance, (if active => adjustments such that there is no trade deficit, i.e. counterfactual in itself)
TB_ctry_prime = copy(TB_ctry)
# TB_ctry_prime .= 0.0


while max_error > tolerance && iteration <= max_iteration

    # Price indices
    P_hat_Z, P_hat_F, π_hat_Z, π_hat_F, cost_hat = create_price_index_hat(w_hat, τ_hat_Z, τ_hat_F, VA_coeff, γ, π_Z, π_F, θ)

    # Counterfactual trade shares
    π_hat_Z = ifelse.(isinf.(π_hat_Z), 0.0, π_hat_Z) # remove Inf
    π_hat_F = ifelse.(isinf.(π_hat_F), 0.0, π_hat_F)

    global π_prime_Z = π_Z .* π_hat_Z # NS×NS
    global π_prime_F = π_F .* π_hat_F # NS×N

    # Labor market clearing
    w_hat_prev = copy(w_hat) # store last wage in case new optimization obtains negative wages

    w_hat, Y_prime, Z_prime, F_prime, ETB_ctry = create_wages_hat(w_hat, vfactor, π_prime_Z, π_prime_F, VA_ctry, TB_ctry_prime, γ, α)

    # update iteration parameters
    error = abs.(w_hat .- w_hat_prev)
    max_error = maximum(error) # update error
    iteration += 1 # update iteration count

    if minimum(w_hat) < 0.0
        w_hat[:] = copy(w_hat_prev)
        println("Outer loop: Iteration $iteration completed with error $max_error (wage negative, rerun with previous estimate)")
    else
        println("Outer loop: Iteration $iteration completed with error $max_error")
    end


end

w_hat

# Country consumer price index and real wage
P0_F_ctry = (ones(S,N) ./ α).^α # S×N
P0_F_ctry = [prod(P0_F_ctry[:,j]) for j in 1:N] # N×1, initial consumer price index

P_F_ctry = (P_hat_F ./ α).^α # S×N, counterfactual consumer price index
P_F_ctry = [prod(P_F_ctry[:,j]) for j in 1:N] # N×1

P_hat_F_ctry = P_F_ctry ./ P0_F_ctry # N×1, changes in country price index

w_hat_real = w_hat ./ P_hat_F_ctry # N×1, real wage changes

# -------------------------------------------------------------------------------------------------------------------------------------

# Gross output
Y_ctry = [sum(Y[i:i+S-1]) for i in 1:S:N*S] # N×1
Y_prime_ctry = [sum(Y_prime[i:i+S-1]) for i in 1:S:N*S] # N×1

Y_hat_ctry = Y_prime_ctry ./ Y_ctry # N×1, nominal

# how to obtain real change since Y is made out of final and intermediate goods?
Y_hat_ctry_real = Y_hat_ctry ./ P_hat_F_ctry