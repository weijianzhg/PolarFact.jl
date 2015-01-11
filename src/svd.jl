#
# Compute Polar Decomposition via SVD
#

immutable SVDAlg <: PolarAlg end

function solve!(alg::SVDAlg, 
                X::Matrix{Float64}, U::Matrix{Float64}, H::Matrix{Float64})

    F = svdfact(X, thin = true)
    PQt = Array(Float64, size(U))

    A_mul_Bt!(PQt, F[:U], F[:V])
    copy!(U, PQt)
    copy!(H, F[:V] * diagm(F[:S]) * F[:Vt])
    
    H = 0.5 * (H + H')
    return Result(U, H)
    
end
