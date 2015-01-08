#
# Compute Polar Decomposition via SVD
#

type SVDAlg
    verbose::Bool

    function SVDAlg(;verbose::Bool=false)
        new(verbose)
    end
end

function solve!(alg::SVDAlg, 
                X::Matrix{Float64}, U::Matrix{Float64}, H::Matrix{Float64})
    # naive implementation
    P, S, Q = svd(X, thin = true)
    U[:] = P * Q'
    H[:] = Q * diagm(S) * Q' 
    H = (H + H') / 2
    return Result(U, H)
    
end
