


using DataFrames, RData, LinearAlgebra, Statistics

years = 2000:2014 # vector with years covered by WIOD Rev. 2016
N = 41 # number of countries 
S = 35 # number of industries
dir = "C:/Users/u0148308/Desktop/raw/" # location of raw data

# --------------- Import data ---------------------------------------------------------------------------------------------------------------------------------

function import_excel(dir::String, year::Integer, N::Integer, S::Integer)

    # load data for specified year and transform into a DataFrame (automatically)
    path = dir * "WIOT" * string(year)[end-1:end] * "_ROW_Apr12.xlsx"
    df = DataFrame(XLSX.readxlsx(path)["WIOT_$year"][:], :auto)
    df = df[7:end,5:end]

    # Intermediate demand
    Z = df[1:N*S, 1:N*S] # NS×NS
    Z = Matrix(convert.(Float64, Z))

    # Final demand
    F = df[1:N*S, N*S+1:end-1] # NS×5*N
    
    # Subset final demand for inventory adjustments
    inventory_columns = 5:5:5*N
    IV = F[:, inventory_columns] # NS×N
    F = F[:, Not(inventory_columns)] # NS×4*N
    
    IV = Matrix(convert.(Float64, IV)) # NS×N
    F = Matrix(convert.(Float64, F)) # NS×4*N


    return Z, F, IV
end

function antras_trans1(Z::Matrix, F::Matrix, IV::Matrix, N::Integer, S::Integer)

    XX = Z
    FF = F
    NN = IV

    NN = [sum(NN[i,:]) for i in 1:N*S] # NS×1
    
    YY = [sum(XX[i,:]) + sum(FF[i,:]) + NN[i] for i in 1:N*S] # NS×1
    YY = ifelse.(abs.(YY) .< 1e-6, 0.0, YY)

    adj = YY.-NN
    adj = ifelse.(adj .< 0.0, 0.0, adj)

    XXi = XX .* (repeat(YY,1,N*S) ./ repeat(adj,1,N*S)) # NS×NS, inventory adjustment
    XXi = ifelse.(isnan.(XXi), 0.0, XXi)
    XXi = ifelse.(isinf.(XXi), 0.0, XXi)

    FFi = FF .* (repeat(YY,1,N*4) ./ repeat(adj,1,N*4)) # NS×N*4, inventory adjustment
    FFi = ifelse.(isnan.(FFi), 0.0, FFi)
    FFi = ifelse.(isinf.(FFi), 0.0, FFi)
    FFdi = [sum(FFi[i,j:j+3]) for i in 1:N*S, j in 1:4:N*4] # NS×N

    YYi = [sum(XXi[i,:]) + sum(FFdi[i,:]) for i in 1:N*S] # NS×1

    VVi = [YYi[i] - sum(XXi[:,i]) for i in 1:N*S] # NS×1

    return XXi, FFdi, YYi, VVi
end


function antras_trans2(XXi::Matrix, FFdi::Matrix, YYi::Vector, VVi::Vector)
    
    gammas = zeros(S+1,N*S)

    for n in 1:N
    for s in 1:S
        inputlist = XXi[:,(n-1)*S+s]
        inputlist = reshape(inputlist, S, N)
        inputlist = [sum(inputlist[i,:]) for i in 1:S]
        inputlist = [inputlist; VVi[(n-1)*S+s]]
        inputratio = inputlist ./ sum(inputlist)
        gammas[:,(n-1)*S+s] = inputratio
    end
    end

    VAcoeff = gammas[S+1,:]
    VAcoeff = ifelse.(isnan.(VAcoeff), 1.0, VAcoeff)
    VAcoeff = reshape(VAcoeff, S, N)

    gammas = gammas[1:S, 1:N*S]
    gammas = ifelse.(isnan.(gammas), 0.0, gammas)

    alphas = zeros(S, N)
    CtyCons = [sum(FFdi[:,i]) for i in 1:N]

    for j in 1:N
        for i in 1:S
            alphas[i, j] = sum(FFdi[i:S:(N-1)*S+i, j]) / CtyCons[j]
        end
    end




    VAn = [sum(VVi[i:i+S-1]) for i in 1:S:N*S] # N×1

    Piagg = [sum(XXi[s:S:(N-1)*S+s,j]) for s in 1:S, j in 1:N*S]
    PiFagg = [sum(FFdi[s:S:(N-1)*S+s,j]) for s in 1:S, j in 1:N]

    Pi = XXi ./ repeat(Piagg, N)
    Pi = ifelse.(isnan.(Pi), 0.0, Pi)

    PiF = FFdi ./ repeat(PiFagg, N)
    PiF = ifelse.(isnan.(PiF), 0.0, PiF)

    M = zeros(N)
    E = zeros(N)

    for j in 1:N
    for s in 1:S
        for i in 1:N
            iiss = (i-1)*S+s
            jjss = (j-1)*S+s
            for r in 1:S
                jjrr = (j-1)*S+r
                iirr = (i-1)*S+r
                M[j] = M[j] + Pi[iiss,jjrr]*gammas[s,jjrr]*YYi[jjrr]
                E[j] = E[j] + Pi[jjss,iirr]*gammas[s,iirr]*YYi[iirr]
            end
            M[j] = M[j] + PiF[iiss,j]*alphas[s,j]*CtyCons[j]
            E[j] = E[j] + PiF[jjss,i]*alphas[s,i]*CtyCons[i]
        end
    end
    end

    Surplusn = E .- M

    return gammas, alphas, Surplusn, VAcoeff, VAn, Pi, PiF, CtyCons
end

Z, F, IV = import_excel(dir, 1995, N, S)
##################
XXi, FFdi, YYi, VVi = antras_trans1(Z, F, IV, N, S)
gammas, alphas, Surplusn, VAcoeff, VAn, Pi, PiF, CtyCons = antras_trans2(XXi, FFdi, YYi, VVi)


Z, F, Y, VA = create_matrices(Z, F, IV, N, S)
π_Z, π_F = create_trade_shares(Z, F, N, S)
γ, α, VA_coeff, TB_ctry, VA_ctry, F_ctry = create_expenditure_shares(Z, F, Y, VA, π_Z, π_F, N, S)
#####################

alphas == α
gammas == γ
PiF == π_F
Pi == π_Z
F == FFdi
F_ctry == CtyCons

VA_ctry.- F_ctry
Surplusn
TB_ctry
Surplusn .- TB_ctry

[sum(α[:,j]) for j in 1:N]
[sum(alphas[:,j]) for j in 1:N]

count(round.(gammas,digits=0) .!= round.(γ,digits=0))

count(abs.(alphas .- α) .> 0.0001)

γ2 = zeros(S+1,N*S)
for n in 1:N
    for s in 1:S
        inputlist = Z[:,(n-1)*S+s]
        inputlist = reshape(inputlist, S, N)
        inputlist = [sum(inputlist[i,:]) for i in 1:S]
        inputlist = [inputlist; VA[(n-1)*S+s]]
        inputratio = inputlist ./ sum(inputlist)
        γ2[:,(n-1)*S+s] = inputratio
    end
end

VA2 = γ2[S+1,:]
VA2 = ifelse.(isnan.(VA2), 1.0, VA2)
VA2 = reshape(VA2, S, N)

γ2 = γ2[1:S, 1:N*S]
γ2 = ifelse.(isnan.(γ2), 0.0, γ2)


count(abs.(gammas .- γ2) .> 5)

count(γ.>1)

findall(γ.>=1)
findall(γ.<0)


findall(VA_coeff.>1)

findall(π_Z.<0)
findall(π_F.<0)

findall(Z.<0)
findall(VA_ctry.<0)


E_A = zeros(N)
M_A = zeros(N)
for j in 1:N
    for s in 1:S
        for i in 1:N
            iiss = (i-1)*S+s
            jjss = (j-1)*S+s
            for r in 1:S
                jjrr = (j-1)*S+r
                iirr = (i-1)*S+r
                M_A[j] += π_Z[iiss,jjrr]*γ[s,jjrr]*Y[jjrr]
                E_A[j] += π_Z[jjss,iirr]*γ[s,iirr]*Y[iirr]
            end
            M_A[j] += π_F[iiss,j]*α[s,j]*F_ctry[j]
            E_A[j] += π_F[jjss,i]*α[s,i]*F_ctry[i]
        end
    end
end

TB_ctry2 = E_A .- M_A # N×1

F_prime_ctry = copy(F_ctry)
Y_prime = copy(Y_prime)

LHS = zeros(N)
RHS = zeros(N)
for j in 1:N
    for s in 1:S
        for i in 1:N
            iiss = (i-1)*S+s
            jjss = (j-1)*S+s
            for r in 1:S
                jjrr = (j-1)*S+r
                iirr = (i-1)*S+r
                LHS[j] += π_Z[iiss,jjrr]*γ[s,jjrr]*Y[jjrr]
                RHS[j] += π_Z[jjss,iirr]*γ[s,iirr]*Y[iirr]
            end
            LHS[j] += π_F[iiss,j]*α[s,j]*F_ctry[j]
            RHS[j] += π_F[jjss,i]*α[s,i]*F_ctry[i]
        end
    end
end

ETB_ctry = RHS .- LHS

findall(abs.(b).>1)



Z_prime = π_Z .* repeat(γ, N) .* repeat(transpose(Y), N*S) # NS×NS
F_prime = π_F .* repeat(α, N) .* repeat(transpose(F_ctry), N*S) # NS×N

E_prime = [sum(Y[i:i+S-1]) - sum(Z_prime[i:i+S-1,i:i+S-1]) - sum(F_prime[i:i+S-1,ceil(Int, i/S)]) for i in 1:S:N*S] # N×1

M_prime_Z = [sum(Z_prime[:,j:j+S-1]) - sum(Z_prime[j:j+S-1,j:j+S-1]) for j in 1:S:N*S]
M_prime_F = [sum(F_prime[:,ceil(Int, j/S)]) - sum(F_prime[j:j+S-1,ceil(Int, j/S)]) for j in 1:S:N*S]
M_prime = M_prime_Z .+ M_prime_F # N×1

ETB_ctry_prime = E_prime .- M_prime

findall(abs.(Z.-Z_prime).>0.1)
findall(abs.(ETB_ctry.-ETB_ctry_prime).>1)

TB_ctry

ETB = ETB_ctry .- TB_ctry

b = ETB./VA_ctry

δ = sign.(b) .* abs.(vfactor.*b)
w_hat = ones(N)
w_hat = w_hat .* (1.0 .+ δ ./ w_hat)


a = TB_ctry./VA_ctry
δ = sign.(a) .* abs.(vfactor.*a)
w_hat = ones(N)
w_hat = w_hat .* (1.0 .+ δ ./ w_hat)

findall(abs.(a).>1)

findall(TB_ctry .> VA_ctry)