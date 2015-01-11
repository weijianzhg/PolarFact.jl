#
# Compute Polar Decomposition via SVD
#

type SVDAlg end

function solve!(alg::SVDAlg, 
                X::Matrix{Float64}, U::Matrix{Float64}, H::Matrix{Float64})

    P, S, Q = svd(X, thin = true)
    PQt = Array(Float64, size(U))
    Ht = Array(Float64, size(U))

    A_mul_Bt!(PQt, P, Q)
    copy!(U, PQt)
    copy!(H, Q * diagm(S) * Q')
    
    copy!(Ht, H')
    for i in length(H)
        H[i] = ( H[i] + Ht[i] )/2
    end
    return Result(U, H)
    
end
