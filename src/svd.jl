#
# Compute Polar Decomposition via SVD
#

immutable SVDAlg <: PolarAlg end

function solve!(alg::SVDAlg, 
                X::Matrix{Float64}, U::Matrix{Float64}, H::Matrix{Float64})

    P, S, Q = svd(X, thin = true)
    PQt = Array(Float64, size(U))

    A_mul_Bt!(PQt, P, Q)
    copy!(U, PQt)
    copy!(H, Q * diagm(S) * Q')
    
    H = 0.5 * (H + H')
    return Result(U, H)
    
end
