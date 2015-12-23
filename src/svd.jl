#
# Compute Polar Decomposition via SVD
#

immutable SVDAlg <: PolarAlg end

function solve!{T}(alg::SVDAlg, 
                X::Matrix{T}, U::Matrix{T}, H::Matrix{T})

    F = svdfact(X, thin = true)
    PQt = Array(T, size(U))

    A_mul_Bt!(PQt, F[:U], F[:V])
    copy!(U, PQt)
    copy!(H, F[:V] * diagm(F[:S]) * F[:Vt])
    H = 0.5 * (H + H')
    
    return SVDResult{T}(U, H)
end
