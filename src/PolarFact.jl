module PolarFact

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
                   tol::T = convert(T, 1.0e-6),
                   verbose::Bool=false)

    # choose algorithm 
    algorithm = 
       alg == :newton ? NewtonAlg(maxiter=maxiter, tol=tol, verbose=verbose) :
       alg == :qdwh ?  QDWHAlg(maxiter=maxiter, tol=tol, verbose=verbose) : 
       alg == :halley ? HalleyAlg(maxiter=maxiter, tol=tol, verbose=verbose) :
       alg == :svd ? SVDAlg() :   
       alg == :schulz ? NewtonSchulzAlg(maxiter=maxiter, tol=tol, verbose=verbose):
       alg == :hybrid ? NewtonHybridAlg(maxiter=maxiter, tol=tol, verbose=verbose) :
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
