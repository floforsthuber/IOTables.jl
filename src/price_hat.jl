# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script to obtain the price index
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

using DataFrames, RData

include("transform_WIOD_2016.jl") # Script with functions to import and transform raw data

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

years = 2000:2014 # vector with years covered by WIOD Rev. 2016
N = 44 # number of countries 
S = 56 # number of industries
dir = "C:/Users/u0148308/Desktop/raw/" # location of raw data

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

Z, F, Y, F_ctry, TB_ctry, VA_ctry, VA_coeff, γ, α, π_Z, π_F = transform_WIOD_2016(dir, 2014)

# further definitions needed (use same names as Antras and Chor (2018) for the time being)

# iteration parameters
vfactor = 0.4
tol = 1e-6
maxit = 500
wfmax = 1e7

# sectoral trade elasticity
cons = 5
θ = cons .* ones(S)

# trade costs
τ_hat_Z = ones(N*S, N*S)
τ_hat_F = ones(N*S, N)

# final country level trade balance
# should take this to be either 0 or the final value in 2011?
TB_new = TB_ctry

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
