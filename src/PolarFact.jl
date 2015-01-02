module PolarFact

include("common.jl") # common functions and types
include("newton.jl") # using Newton's method
# include("newton_hybrid.jl")
include("svd.jl")
# include("")

export polarfact

function polarfact{T}(A::Matrix{T}; 
                      alg::Symbol=:newton, 
                      maxiter::Integer=:100,
                      tol::Real=1.0e-6,
                      verbose::Bool=false)

    # choose algorithm 
    algorithm = 
    alg == :newton ? NewtonAlg(maxiter=maxiter, tol=tol, verbose=verbose) :
    error("Invalid algorithm.")

    # Initialization: if m > n, do QR factorization 
    m, n = size(A)
    if m > n
        qrfact!(A)
    end

    U = Array(Float64, size(A))
    H = Array(Float64, size(A))
    # solve for polar factors
    solve!(algorithm, A, U, H)
end

end # module
