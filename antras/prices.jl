

S = 56
N = 44

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