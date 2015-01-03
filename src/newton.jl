#
# Newton's method
#
# Reference:
# Nicholas J. Higham, Computing the Polar Decomposition ---with Applications,
# SIAM J. Sci. Statist. Comput. Vol. 7, Num 4 (1986) pp. 1160-1174.
#
   
type NewtonAlg
    maxiter::Int        # maximum number of iterations.
    scale::Bool         # whether to scale Newton iteration.
    verbose::Bool       # whether to show procedural information
    tol::Float64        # tolerance for convergence
    scale_tol::Float64  # tolerance for acceleration scaling 

    function NewtonAlg(;maxiter::Integer=100,
                       verbose::Bool=false,
                       scale::Bool=true,
                       tol::Real=1.0e-6,
                       scale_tol::Real=1.0e-2)
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        scale_tol > 0 || error("scale_tol must be positive.")
        
        new(int(maxiter),
            scale,
            verbose,
            float64(tol),
            float64(scale_tol))
    end
end


function solve!(alg::NewtonAlg, 
                X::Matrix{Float64}, U::Matrix{Float64}, H::Matrix{Float64})
    common_iter!(NewtonUpdater(alg.scale, alg.scale_tol), X, U, H, alg.maxiter, alg.verbose, alg.tol)
end

type NewtonUpdater <: PolarUpdater 
    scale::Bool
    scale_tol::Float64
end

function update_U!(upd::NewtonUpdater, U::Matrix{Float64})
    scale = upd.scale
    Uinv = Array(Float64, size(U))
    Uinvt = Array(Float64, size(U))
    Uinv = inv(U)
    
    # 1, Inf-norm scaling 
    if scale
        g = (norm(Uinv,1) * norm(Uinv, Inf) / (norm(U,1) * norm(U, Inf)) )^(1/4)
    else
        g = 1
    end
    transpose!(Uinvt, Uinv)
    for i = 1:length(U)
        U[i] = (g * U[i] + Uinvt[i] / g)/2
    end

end


