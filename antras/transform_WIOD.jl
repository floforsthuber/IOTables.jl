


using DataFrames, RData, LinearAlgebra, Statistics

years = 2000:2014 # vector with years covered by WIOD Rev. 2016
N = 44 # number of countries 
S = 56 # number of industries
dir = "C:/Users/u0148308/Desktop/raw/" # location of raw data

# --------------- Import data ---------------------------------------------------------------------------------------------------------------------------------

function import_R(dir::String, year::Integer)

    # load data for specified year and transform into a DataFrame (automatically)
    path = dir * "WIOT" * string(year) * "_October16_ROW.RData"
    df = RData.load(path)["wiot"] # 2472×2690 DataFrame
    transform!(df, [:Year, :RNr] .=> ByRow(Int64) .=> [:Year, :RNr], renamecols=false) # Years, row-industry-identifier as type: Int64

    return df
end


function antras_trans1(df::DataFrame, N::Integer, S::Integer)

    XX = df[1:N*S, 6:N*S+5] # NS×NS
    XX = Matrix(convert.(Float64, XX))

    FF = df[1:N*S, N*S+6:end-1] # NS×5*N
    
    inventory_columns = 5:5:5*N
    NN = FF[:, inventory_columns] # NS×N
    FF = FF[:, Not(inventory_columns)] # NS×4*N
    
    NN = Matrix(convert.(Float64, NN)) # NS×N
    FF = Matrix(convert.(Float64, FF)) # NS×4*N

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
    FFdi = [sum(FF[i,j:j+3]) for i in 1:N*S, j in 1:4:N*4] # NS×N

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


##################
df = import_R(dir, 2014)
XXi, FFdi, YYi, VVi = antras_trans1(df, N, S)
gammas, alphas, Surplusn, VAcoeff, VAn, Pi, PiF, CtyCons = antras_trans2(XXi, FFdi, YYi, VVi)


Z, F, Y, VA = create_matrices(df, N, S)
π_Z, π_F = create_trade_shares(Z, F, N, S)
γ, α, VA_coeff, TB_ctry, VA_ctry, F_ctry = create_expenditure_shares(Z, F, Y, VA, π_Z, π_F, N, S)
#####################

count(round.(gammas,digits=2) .!= round.(γ,digits=2))

count(round.(alphas,digits=1) .!= round.(α,digits=1))
