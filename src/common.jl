# common types or functions

# the result type

immutable Result{T}
    U::Matrix{T}
    H::Matrix{T}
    niters::Int
    converged::Bool
end

# common algorithm skeleton for iterative updating methods

function common_alg!(updater::PolarUpdater,
                     X::Matrix{Float64}, 
                     U::Matrix{Float64}, 
                     H::Matrix{Float64},
                     maxiter::Int,
                     verbose::Bool,
                     tol::Float64)
    
    preU = Array(Float64, size(U))
    
    converged = false
    t = 0
    while !converged && t < maxiter
        t +=1
        copy!(preU, U)
        
        update_UH!(updater, X, U)
        
        # determine convergence
        mu = norm()

    end
end
