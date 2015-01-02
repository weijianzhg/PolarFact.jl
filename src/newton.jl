#
# Newton's method
#
# Reference:
# Nicholas J. Higham, Computing the Polar Decomposition ---with Applications,
# SIAM J. Sci. Statist. Comput. Vol. 7, Num 4 (1986) pp. 1160-1174.
#
   
type NewtonAlg
    maxiter::Int      # maximum number of iterations.
    verbose::Bool     # whether to show procedural information
    tol::Float64      # tolerance

    function NewtonAlg(;maxiter::Integer=100,
                       verbose::Bool=false,
                       tol::Real=1.0e-6)
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive")
        
        new(int(maxiter), 
            verbose,
            float64(tol))
    end
end


function solve!(alg::NewtonAlg, 
                X::Matrix{Float64}, U::Matrix{Float64}, H::Matrix{Float64})
    common_iter!(NewtonUpdater(), X, U, H, alg.maxiter, alg.verbose, alg.tol)
end

immutable NewtonUpdater <: PolarUpdater end

function update_U!(upd::NewtonUpdater, U::Matrix{Float64})
    invU = Array(Float64, size(U))
    invU = inv(U)
    for i = 1:length(U)
        U[i] = (U[i] + invU[i])/2
    end
end
