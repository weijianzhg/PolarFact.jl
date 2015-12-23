module PolarFact

using Compat
import Base: transpose!

export polarfact


include("common.jl")     # common functions and types
include("newton.jl")     # using Newton's method
include("svd.jl")        # using the SVD method
include("halley.jl")     # using Halley's method and the QDWH method
include("hybrid.jl")     # using hybrid method

function polarfact{T}(A::Matrix{T}; 
                   alg::Symbol=:newton, 
                   maxiter::Integer=:100,
                   tol::Real = cbrt(eps(T)),
                   verbose::Bool=false)

    # choose algorithm 
    algorithm = 
       alg == :newton ? NewtonAlg{T}(maxiter=maxiter, tol=tol, verbose=verbose) :
       alg == :qdwh ?  QDWHAlg{T}(maxiter=maxiter, tol=tol, verbose=verbose) : 
       alg == :halley ? HalleyAlg{T}(maxiter=maxiter, tol=tol, verbose=verbose) :
       alg == :svd ? SVDAlg() :   
       alg == :schulz ? NewtonSchulzAlg{T}(maxiter=maxiter, tol=tol, verbose=verbose):
       alg == :hybrid ? NewtonHybridAlg{T}(maxiter=maxiter, tol=tol, verbose=verbose) :
       error("Invalid algorithm.")

    # Initialization: if m > n, do QR factorization 
    m, n = size(A)
    if m > n
        A = qrfact(A)[:R]
    elseif m < n
        error("The row dimension of the input matrix must be 
              greater or equal to column dimension.")
    end

    U = Array(T, size(A))
    H = Array(T, size(A))
    # solve for polar factors
    solve!(algorithm, A, U, H)
end

end # module
