# 
# Halley's method
# 
# Reference:
# Y. Nakatsukasa, Z. Bai and F. Gygi, Optimizing Halley's iteration 
# for computing the matrix polar decomposition, SIAM, J. Mat. Anal. 
# Vol. 31, Num 5 (2010) pp. 2700-2720 
#

type HalleyAlg <: PolarAlg
    maxiter::Int
    verbose::Bool
    tol::Float64
    
    function HalleyAlg(;maxiter::Integer=100,
                       verbose::Bool=false,
                       tol::Real=1.0e-6)
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        
        new(int(maxiter),
            verbose,
            float64(tol))
    end
end

function solve!(alg::HalleyAlg,
                X::Matrix{Float64}, U::Matrix{Float64}, H::Matrix{Float64})
    common_iter!(HalleyUpdater(), X, U, H, alg.maxiter, alg.verbose, alg.tol)
end

immutable HalleyUpdater <: PolarUpdater end


function update_U!(upd::HalleyUpdater, U::Matrix{Float64})   
    UtU = Array(Float64, size(U))
    At_mul_B!(UtU, U, U)
    copy!(U, U * (3*I + UtU)* inv(I + 3*UtU))
end

#
# QR-based Dynamically Weighted Halley (QDWH) algorithm
#

type QDWHAlg <: PolarAlg
    maxiter::Int
    verbose::Bool
    piv::Symbol  # type of pivoting
    tol::Float64
    
    function QDWHAlg(;maxiter::Integer=100,
                     verbose::Bool=false,
                     piv::Symbol= :rc, # column pivoting and row sorting
                     tol::Real=1.0e-6)
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        
        new(int(maxiter),
            verbose,
            piv,
            float64(tol))
    end
end


function solve!(alg::QDWHAlg,
                X::Matrix{Float64}, U::Matrix{Float64}, H::Matrix{Float64})
    common_iter!(QDWHUpdater(alg.piv), X, U, H, alg.maxiter, alg.verbose, alg.tol)
end

type QDWHUpdater <: PolarUpdater
    piv::Symbol
    alpha::Float64
    L::Float64
end


function update_U!(upd::QDWHUpdater, U::Matrix{Float64})   
    piv = upd.piv
    alpha = upd.alpha
    L = upd.L

end

