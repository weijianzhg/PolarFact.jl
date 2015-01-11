# common types or functions

# the result type

immutable Result
    U::Matrix{Float64}
    H::Matrix{Float64}
    niters::Union(Int, Nothing)
    converged::Union(Bool, Nothing)
    
    function Result(U::Matrix{Float64}, H::Matrix{Float64}, 
                    niter::Union(Int,Nothing)=nothing, 
                    converged::Union(Bool, Nothing)=nothing)
        size(U, 2) == size(H, 1) || 
               throw(DimensionMismatch("Inner dimension of U and H mismatch."))
        new(U, H, niter, converged)
    end
end

# niters and converged are nothing for SVD method

# the objective type

immutable Objective
    absolute::Float64  # absolute error
    relative::Float64  # relative error

    function Objective(absolute::Real, relative::Real)
        absolute >= 0 || error("absolute must be non-negative.")
        relative >= 0 || error("relative must be non-negative.")

        new(float64(absolute), float64(relative))
    end
end


function evaluate_objv(preU::Matrix{Float64}, U::Matrix{Float64})
    rel = norm(preU - U, Inf) / norm(preU, Inf)
    abs = norm(I - U'*U, Inf)
    return Objective(rel, abs)
end


abstract PolarUpdater

abstract PolarAlg

# common algorithm skeleton for iterative updating methods

function common_iter!(updater::PolarUpdater,
                      X::Matrix{Float64}, 
                      U::Matrix{Float64}, 
                      H::Matrix{Float64},
                      maxiter::Int,
                      verbose::Bool,
                      tol::Float64)
    
    preU = Array(Float64, size(X))
    copy!(U, X)
    converged = false
    t = 0
    if verbose
        @printf("%-5s    %-13s    %-13s\n", "Iter.", "Rel. err.",  "Obj.")
    end

    while !converged && t < maxiter
        t += 1
        copy!(preU, U)
        update_U!(updater, U)
        
        # determine convergence
        diff = vecnorm(preU - U)
        if diff < tol
            converged = true
        end
        
        # display infomation
        if verbose
            objv = evaluate_objv(preU, U)
            @printf("%5d    %13.6e    %13.6e\n", 
                    t, objv.absolute, objv.relative)
        end
    end

    # compute H
    A_mul_B!(H, U', X)
    H = 0.5 * (H + H')
    return Result(U, H, t, converged)
end

# Scaling iterative algorithm
function common_iter_scal!(updater::PolarUpdater,
                            X::Matrix{Float64}, 
                            U::Matrix{Float64}, 
                            H::Matrix{Float64},
                            maxiter::Int,
                            verbose::Bool,
                            tol::Float64)
    
    preU = Array(Float64, size(X))
    copy!(U, X)
    converged = false
    t = 0
    if verbose
        @printf("%-5s    %-13s    %-13s\n", "Iter.", "Rel. err.",  "Obj.")
    end

    while !converged && t < maxiter
        t += 1
        copy!(preU, U)
        update_U!(updater, U)
        
        # determine convergence
        diff = vecnorm(preU - U)
        if diff < tol
            converged = true
        end
        
        # determine scaling
        reldiff = diff/vecnorm(U) # relative error
        if updater.scale && (reldiff < updater.scale_tol)
            updater.scale = false
        end

        # display infomation
        if verbose
            objv = evaluate_objv(preU, U)
            @printf("%5d    %13.6e    %13.6e\n", 
                    t, objv.absolute, objv.relative)
        end
    end

    # compute H
    A_mul_B!(H, U', X)
    H = 0.5 * (H + H')
    return Result(U, H, t, converged)
end



# Hybrid iteration algorithm 
function common_iter_hybr!(updater1::PolarUpdater,
                           updater2::PolarUpdater,
                           X::Matrix{Float64}, 
                           U::Matrix{Float64}, 
                           H::Matrix{Float64},
                           maxiter::Int,
                           verbose::Bool,
                           tol::Float64,
                           theta::Float64)

    preU = Array(Float64, size(X))
    copy!(U, X)
    converged = false
    switched = false
    t = 0
    if verbose
        @printf("%-5s    %-13s    %-13s\n", "Iter.", "Rel. err.",  "Obj.")
    end

    while !converged && t < maxiter
        t += 1
        copy!(preU, U)

        if switched
            update_U!(updater2, U)        
        else
            obj = norm(I - U'*U, 1)
            if obj > theta # theta is the switch parameter
                update_U!(updater1, U)
            else
                switched = true
                update_U!(updater2, U)
            end
        end

        # determine convergence
        diff = vecnorm(preU - U)
        if diff < tol
            converged = true
        end
        
        # display infomation
        if verbose
            objv = evaluate_objv(preU, U)
            @printf("%5d    %13.6e    %13.6e\n", 
                    t, objv.absolute, objv.relative)
        end
    end

    # compute H
    A_mul_B!(H, U', X)
    H = 0.5 * (H + H')
    return Result(U, H, t, converged)

end
