#
# Newton's method
#
# Reference:
# Nicholas J. Higham, Computing the Polar Decomposition ---with Applications,
# SIAM J. Sci. Statist. Comput. Vol. 7, Num 4 (1986) pp. 1160-1174.
#
   
type NewtonAlg{T} <: PolarAlg
    maxiter::Int        # maximum number of iterations.
    scale::Bool         # whether to scale Newton iteration.
    verbose::Bool       # whether to show procedural information
    tol::T        # tolerance for convergence
    scale_tol::T  # tolerance for acceleration scaling 

    function NewtonAlg(;maxiter::Integer=100,
                       verbose::Bool=false,
                       scale::Bool=true,
                       tol::Real=cbrt(eps(T)),
                       scale_tol::Real=eps(T)^(1/4))
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        scale_tol > 0 || error("scale_tol must be positive.")
        
        new(maxiter,
            scale,
            verbose,
            tol,
            scale_tol)
    end
end


function solve!{T}(alg::NewtonAlg{T}, 
                   X::Matrix{T}, U::Matrix{T}, H::Matrix{T})
    common_iter_scal!(NewtonUpdater(alg.scale, alg.scale_tol), X, U, H, alg.maxiter, alg.verbose, alg.tol)
end

type NewtonUpdater{T} <: PolarUpdater 
    scale::Bool
    scale_tol::T
end

function update_U!{T}(upd::NewtonUpdater, U::Matrix{T})
    scale = upd.scale
    Uinv = Array(T, size(U))
    Uinvt = Array(T, size(U))
    copy!(Uinv, inv(U))
    
    # 1, Inf-norm scaling 
    if scale
        g = (norm(Uinv,1) * norm(Uinv, Inf) / (norm(U,1) * norm(U, Inf)) )^(1/4)
    else
        g = one(T)
    end
    transpose!(Uinvt, Uinv)
    copy!(U, (g * U + Uinvt /g) / convert(T, 2))

end



# Newton Schulz algrithm

type NewtonSchulzAlg{T} <: PolarAlg  
    maxiter::Int     # maximum number of iterations.
    verbose::Bool    # whether to show procedural information
    tol::T     # convergence tolerance. 

    function NewtonSchulzAlg(;maxiter::Integer=100,
                             verbose::Bool=false,
                             tol::Real=cbrt(eps(T)))
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")

        new(int(maxiter),
            verbose,
            tol)
    end
end

function solve!{T}(alg::NewtonSchulzAlg, 
                   X::Matrix{T}, U::Matrix{T}, H::Matrix{T})
    
    # Newton Schulz converge quadratically if norm(X) < sqrt(3)
    norm(X) < convert(T, sqrt(3)) || throw(ArgumentError("The norm of the input matrix must be smaller than sqrt(3)."))

    common_iter!(NewtonSchulzUpdater(), X, U, H, alg.maxiter, alg.verbose, alg.tol)
end

immutable NewtonSchulzUpdater <: PolarUpdater end

function update_U!{T}(upd::NewtonSchulzUpdater, U::Matrix{T})
    UtU = Array(T, size(U))
    At_mul_B!(UtU, U, U)
    copy!(U, 0.5*U*(3*I - UtU))
end
