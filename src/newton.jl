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
    theta::Float64    # parameter to change to Newton Schulz iteration

    function NewtonAlg(maxiter::Integer=100,
                       verbose::Bool=false,
                       tol::Real=1.0e-6,
                       theta::Real=0.6)
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive")
        theta >= 0 || error("theta must be non-negative")
        
        new(int(maxiter), 
            verbose,
            float64(tol),
            float64(theta))
    end
end

function solve!(alg::NewtonAlg, 
                X::Matrix{Float64}, U::Matrix{Float64}, H::Matrix{Float64})

end
