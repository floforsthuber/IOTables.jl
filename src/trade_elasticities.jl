# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Baseline model
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData,  XLSX, LinearAlgebra, Statistics

include("transform_data2.jl") # Script with functions to import and transform raw data
include("price_hat.jl") # Script with function to obtain the price index
include("wage_hat.jl") # Script with function to obtain the wages and gross output
include("tariffs_function.jl") # Script with functions to create τ_hat_Z, τ_hat_F from tariff data
include("head_ries_index.jl") # Script with functions to create bilateral Head-Ries index (symmetric bilateral trade costs)
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

years = 2000:2014 # vector with years covered by WIOD Rev. 2016
N = 41 # number of countries 
S = 35 # number of industries
dir = "C:/Users/u0148308/Desktop/raw/" # location of raw data

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_data(dir, "2013", 1995, N, S)

# still need sectoral tariff data (can be asymmetric)
# use Head-Ries Index for the time being (same format, i.e. no internal tariffs)
# notice I assume trade elasticity θ in order to compute them, whilst actual goal is to estimate θ!
θ = 5
τ_HR_Z, τ_HR_F = head_ries_index(Z, F, θ, N, S)


# concord to country-industry exports to country (X_in^j)
X = [sum(Z[i,j:j+S-1]) for i in 1:N*S, j in 1:S:N*S] # NS×N, notice F is missing => final demand specific trade elasticities
τ_X = [mean(τ_HR_Z[i,j:j+S-1]) for i in 1:N*S, j in 1:S:N*S] # NS×N
τ_F = τ_HR_F # NS×N






function create_reg_data(X, F, τ_X, τ_F)

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

    combinations = Int64(sum([n*(n+1)/2 for n in 1:N-2])) # number of possible combinations from formula
    ([length(LHS_Z) length(RHS_Z) length(LHS_F) length(RHS_F)] .== combinations) == [1 1 1 1]

    # create series in logs
    lhs_Z_temp = log.(LHS_Z)
    rhs_Z_temp = log.(RHS_Z)

    lhs_F_temp = log.(LHS_F)
    rhs_F_temp = log.(RHS_F)


    return lhs_Z_temp, rhs_Z_temp, lhs_F_temp, rhs_F_temp
end

lhs_Z = Array{Any, 1}[]
rhs_Z = Array{Any, 1}[]

lhs_F = Array{Any, 1}[]
rhs_F = Array{Any, 1}[]


for j in 1:S
    X_temp = X[j:S:(N-1)*S+j, :]
    τ_X_temp = τ_X[j:S:(N-1)*S+j, :]
    F_temp = F[j:S:(N-1)*S+j, :]
    τ_F_temp = τ_F[j:S:(N-1)*S+j, :]

    lhs_Z_temp, rhs_Z_temp, lhs_F_temp, rhs_F_temp = create_reg_data(X_temp, F_temp, τ_X_temp, τ_F_temp)
    push!(lhs_Z, lhs_Z_temp)
    push!(rhs_Z, rhs_Z_temp)
    push!(lhs_F, lhs_F_temp)
    push!(rhs_F, rhs_F_temp)
end

# zero trade causes inf => create NaN instead which we let regression filter out later
lhs_Z = ifelse.(isinf.(lhs_Z), NaN, lhs_Z)
rhs_Z = ifelse.(isinf.(rhs_Z), NaN, rhs_Z)

lhs_F = ifelse.(isinf.(lhs_F), NaN, lhs_F)
rhs_F = ifelse.(isinf.(rhs_F), NaN, rhs_F)
