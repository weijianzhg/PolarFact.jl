#
# Newton and Newton Schulz Hybrid algoritm
#
# Reference:
# Nicholas J. Higham and Robert S. Schreiber, Fast Polar Decomposition
# of an arbitrary matrix, SIAM, J. Sci. Statist. Comput. Vol. 11, No. 4
# (1990) pp. 648-655.
#

type NewtonHybridAlg <: PolarAlg
    maxiter::Int    # maximum number of iterations.
    verbose::Bool   # whether to show procedural information. 
    tol::Float64    # convergence tolerance
    theta::Float64  # switch parameter
    
    function NewtonHybridAlg(;maxiter::Integer=100,
                             verbose::Bool=false,
                             tol::Real=1.0e-6,
                             theta::Real=0.6)
        maxiter > 1 || error("maxiter must  be greater than 1.")
        tol > 0 || error("tol must be positive.")
        theta > 0 || error("theta must be positive.")
        
        new(int(maxiter), 
            verbose,
            float64(tol),
            float64(theta))
    end 
end
    
function solve!(alg::NewtonHybridAlg, 
                X::Matrix{Float64}, U::Matrix{Float64}, H::Matrix{Float64})
    common_iter_hybr!(NewtonUpdater(true, 1.0e-2), NewtonSchulzUpdater(), X, U, H, 
                      alg.maxiter, alg.verbose, alg.tol, alg.theta)
end

