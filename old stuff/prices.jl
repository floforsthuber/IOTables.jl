

S = 35
N = 41

Theta = ones(S)
tauhat = rand(N*S, N*S)
tauFhat = rand(N*S, N)

wagehat = rand(N)
iter = 1


wagehat0 = wagehat
PrsjHat0 = rand(S, N*S)
PrFjHat0 = rand(S, N)

tol = 1e-3
maxit = 100
pfmax = 1
iterph = 1

cost = ones(S, N)

for s in 1:S
    for j in 1:N
            cost[s,j] = wagehat[j]^VAcoeff[s,j]
        for r in 1:S
            cost[s,j] = cost[s,j]*PrsjHat0[r,(j-1)*S+s]^gammas[r,(j-1)*S+s]
        end
    end
end

PrsijHat = zeros(N*S, N*S)

for i in 1:N
    for r in 1:S
        for j in 1:N
            for s in 1:S
                iirr = (i-1)*S+r
                jjss = (j-1)*S+s
                PrsijHat[iirr,jjss] = (tauhat[iirr,jjss]*cost[r,i])^(-Theta[r])
            end
        end
    end
end

PrsijHat



cost_Z = Pi .* PrsijHat
PrsjHat = [sum(cost_Z[i:S:(N-1)*S+i, j])^(-1/θ[i,ceil(Int,j/S)]) for i in 1:S, j in 1:N*S]


PrFijHat = zeros(N*S,N)
for i in 1:N
    for r in 1:S
        for j in 1:N
            iirr = (i-1)*S+r
            PrFijHat[iirr,j] = (tauFhat[iirr,j]*cost[r,i])^(-Theta[r])
        end
    end
end

PrFijHat

############### mine
cost_hat_w = [wagehat[ceil(Int,i/S)]^VA_coeff[i] for i in 1:N*S] # NS×1, wages
cost_hat_p = PrsjHat0.^γ
cost_hat_p = [prod(cost_hat_p[:,i]) for i in 1:N*S] # NS×1, price indices (prod is the same as sum just for multiplication)
cost_hat = cost_hat_w .* cost_hat_p # NS×1
cost_hat2 = reshape(cost_hat, S, N)


θ = ones(N*S)
cost_hat_Z = zeros(N*S, N*S)
θ = reshape(θ,S,N)

for i in 1:N
    for r in 1:S
        for j in 1:N
            for s in 1:S
                iirr = (i-1)*S+r
                jjss = (j-1)*S+s
                cost_hat_Z[iirr,jjss] = (cost_hat2[r,i]*tauhat[iirr,jjss])^(-θ[r,i])
            end
        end
    end
end

count(abs.(cost_hat_Z .- PrsijHat) .> 1e-5)

cost_Z = π_Z .* cost_hat_Z
P_hat_Z = [sum(cost_Z[i:S:(N-1)*S+i, j])^(-1/θ[i,ceil(Int,j/S)]) for i in 1:S, j in 1:N*S]

count(abs.(PrsjHat .- P_hat_Z) .> 1e-15)
count(isinf.(P_hat_Z))

cost_hat_F = zeros(N*S,N)
for i in 1:N
    for r in 1:S
        for j in 1:N
            iirr = (i-1)*S+r
            cost_hat_F[iirr,j] = (cost_hat2[r,i]*tauFhat[iirr,j])^(-θ[r,i])
        end
    end
end
cost_hat_F

count(abs.(PrFijHat .- cost_hat_F) .> 1e-10)

cost_F = π_F .* cost_hat_F

P_hat_F = [sum(cost_F[i:S:(N-1)*S+i, j])^(-1/θ[i,j]) for i in 1:S, j in 1:N] # S×N


θ = vec(θ)
θ = repeat(θ, 1, N*S)
π_hat_Z = (cost_hat_Z ./ repeat(P_hat_Z, N)) .^ (.-θ) # NS×NS



# -----------------------------------------------------------------------------------------------------------------------------
# wages




# Antras
# -----------------


CtyConsp = wagehat .* VAn .- Surplusn
FGS = π_prime_F .* repeat(alphas, N) .* repeat(CtyConsp', N*S)
FGS = [sum(FGS[i,:]) for i in 1:N*S]
IGS = π_prime_Z .* repeat(gammas, N)
Youtput = inv(I - IGS)*FGS

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
                LHS[j] += π_prime_Z[iiss,jjrr]*gammas[s,jjrr]*Youtput[jjrr]
                RHS[j] += π_prime_Z[jjss,iirr]*gammas[s,iirr]*Youtput[iirr]
            end
            LHS[j] += π_prime_F[iiss,j]*alphas[s,j]*CtyConsp[j]
            RHS[j] += π_prime_F[jjss,i]*alphas[s,i]*CtyConsp[i]
        end
    end
end

ExcessSurplus = RHS .- LHS .- Surplusn # N×1, excess trade balance

vfactor = 0.4
NES = ExcessSurplus ./ VAn
delta = sign.(NES) .* abs.(vfactor .* NES)
wagehat_new = wagehat .* (1 .+ delta ./ wagehat)

# MINE ---------------------

F_prime_ctry = wagehat .* VA_ctry .- TB_ctry # N×1, counterfactual country final goods consumption
# Goods market clearing from equation (35)
total_sales_F = π_prime_F .* repeat(α,N) .* repeat(transpose(F_prime_ctry),N*S) # NS×N
total_sales_F = [sum(total_sales_F[i,:]) for i in 1:N*S] # N×1

A_prime = π_prime_Z .* repeat(γ, N) # NS×NS, intermediate input coefficient matrix

Y_prime = inv(I - A_prime) * total_sales_F #  NS×1


count(abs.(Surplusn .- TB_ctry) .> 0.1)

count(abs.(VAn .- VA_ctry) .> 0.1)

count(abs.(Youtput .- Y_prime) .> 0.1)


Z_prime = π_prime_Z .* repeat(γ, N) .* repeat(transpose(Y_prime), N*S) # NS×NS
F_prime = π_prime_F .* repeat(α, N) .* repeat(transpose(F_prime_ctry), N*S) # NS×N

E_prime = [sum(Y_prime[i:i+S-1]) - sum(Z_prime[i:i+S-1,i:i+S-1]) - sum(F_prime[i:i+S-1,ceil(Int, i/S)]) for i in 1:S:N*S] # N×1

M_prime_Z = [sum(Z_prime[:,j:j+S-1]) - sum(Z_prime[j:j+S-1,j:j+S-1]) for j in 1:S:N*S]
M_prime_F = [sum(F_prime[:,ceil(Int, j/S)]) - sum(F_prime[j:j+S-1,ceil(Int, j/S)]) for j in 1:S:N*S]
M_prime = M_prime_Z .+ M_prime_F # N×1

ETB_ctry = E_prime .- M_prime .- TB_ctry # N×1

a = ExcessSurplus .- ETB_ctry

norm_ETB_ctry =  ExcessSurplus ./ VA_ctry # N×1, normalized excess trade balance
δ = sign.(norm_ETB_ctry) .* abs.(vfactor.*norm_ETB_ctry) # i.e. if norm_ETB_ctry > 0 => too much exports in model => wages should increase (sign fⁿ gives ± of surplus)
w_hat = wagehat .* (1.0 .+ δ ./ wagehat) # N×1, increase/decrease wages for countries with an excess surplus/deficit





F_prime_ctry = w_hat .* VA_ctry .- TB_ctry # N×1, counterfactual country final goods consumption
# Goods market clearing from equation (35)
total_sales_F = π_prime_F .* repeat(α,N) .* repeat(transpose(F_prime_ctry),N*S) # NS×N
total_sales_F = [sum(total_sales_F[i,:]) for i in 1:N*S] # NS×1
A_prime = π_prime_Z .* repeat(γ, N) # NS×NS, intermediate input coefficient matrix
Y_prime = inv(I - A_prime) * total_sales_F #  NS×1

ts = π_F .* repeat(α,N) .* repeat(transpose(F_ctry),N*S)
ts = [sum(ts[i,:]) for i in 1:N*S] # NS×1
A = π_Z .* repeat(γ, N)
Y2 = inv(I - A) * ts

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

# adjust wages to excess trade balance
norm_ETB_ctry =  ETB_ctry ./ VA_ctry # N×1, normalized excess trade balance
δ = sign.(norm_ETB_ctry) .* abs.(vfactor.*norm_ETB_ctry) # i.e. if norm_ETB_ctry > 0 => too much exports in model => wages should increase (sign fⁿ gives ± of surplus)
w_hat3 = w_hat .* (1.0 .+ δ ./ w_hat) # N×1, increase/decrease wages for countries with an excess surplus/deficit






w_hat = rand(N)
P0_hat_Z = rand(S, N*S)

cost_hat_w = [w_hat[ceil(Int,i/S)]^VA_coeff[i] for i in 1:N*S] # NS×1, wages
cost_hat_p = P0_hat_Z.^γ
cost_hat_p = [prod(cost_hat_p[:,i]) for i in 1:N*S] # NS×1, price indices (prod is the same as sum just for multiplication)
cost_hat = cost_hat_w .* cost_hat_p # NS×1
cost_hat = reshape(cost_hat, S, N)

cost_Z = π_Z .* cost_hat_Z # NS×NS, origin country-industry destination country-industry price index (inside of summation)

# sum over origin countries => price index of country-industry composite good in industry
P_hat_Z = [sum(cost_Z[i:S:(N-1)*S+i, j])^(-1/θ[i,ceil(Int,j/S)]) for i in 1:S, j in 1:N*S]

######################

cost = ones(S,N)
VAcoeff = reshape(VA_coeff, S, N)
for s in 1:S
    for j in 1:N
        cost[s,j] = w_hat_prev[j]^VAcoeff[s,j]
        for r in 1:S
            cost[s,j] = cost[s,j]*P_hat_Z[r,(j-1)*S+s]^γ[r,(j-1)*S+s]
        end
    end
end

PrsijHat = zeros(N*S, N*S)
θ = reshape(θ, S, N)

for i in 1:N
    for r in 1:S
        for j in 1:N
            for s in 1:S
                iirr = (i-1)*S+r
                jjss = (j-1)*S+s
                PrsijHat[iirr,jjss] = (τ_hat_Z[iirr,jjss]*cost[r,i])^(-θ[r,i])
            end
        end
    end
end

PrsijHat


PrsjHat = zeros(S,N*S)

for r in 1:S
    for j in 1:N
        for s in 1:S
            jjss = (j-1)*S+s
            PrsjHat[r,jjss] = (sum(π_Z[r:S:(N-1)*S+r,jjss] .* PrsijHat[r:S:(N-1)*S+r,jjss]))^(-1/θ[r,j])
        end
    end
end
PrsjHat

PiHat = PrsijHat ./ (repmat(PrsjHat,N,1))








w_hat = ones(N)


F_prime_ctry = w_hat .* VA_ctry .- TB_ctry # N×1, counterfactual country final goods consumption

# Goods market clearing from equation (35)
total_sales_F = π_F .* repeat(α,N) .* repeat(transpose(F_prime_ctry),N*S) # NS×N
total_sales_F = [sum(total_sales_F[i,:]) for i in 1:N*S] # N×1

A_prime = π_Z .* repeat(γ, N) # NS×NS, intermediate input coefficient matrix

Y_prime = inv(I - A_prime) * total_sales_F #  NS×1
Y_prime = ifelse.(Y_prime .< 0.0, 0.0, Y_prime) # Antras and Chor (2018)


F_prime_ctry ≈ F_ctry

count(abs.(F_prime_ctry.-F_ctry).>0.7)


Z_prime = π_Z .* repeat(γ, N) .* repeat(transpose(Y), N*S) # NS×NS
F_prime = π_F .* repeat(α, N) .* repeat(transpose(F_ctry), N*S) # NS×N
count(abs.(Z_prime.-Z).>1)
count(abs.(F_prime.-F).>1)