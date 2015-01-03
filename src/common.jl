# common types or functions

# the result type

immutable Result
    U::Matrix{Float64}
    H::Matrix{Float64}
    niters::Int
    converged::Bool
end


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
    H = U'*X
    H = 0.5 * (H + H')
    return Result(U, H, t, converged)
end

# common hybrid iteration skeleton 
function common_hybrid!(updater1::PolarUpdater,
                        updater2::PolarUpdater,
                        X::Matrix{Float64}, 
                        U::Matrix{Float64}, 
                        H::Matrix{Float64},
                        maxiter::Int,
                        verbose::Bool,
                        tol::Float64)

end
