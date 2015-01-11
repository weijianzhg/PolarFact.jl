#
# Newton and Newton Schulz Hybrid algoritm
#
# Reference:
# Nicholas J. Higham and Robert S. Schreiber, Fast Polar Decomposition
# of an arbitrary matrix, SIAM, J. Sci. Statist. Comput. Vol. 11, No. 4
# (1990) pp. 648-655.
#

type NewtonHybrid <: PolarAlg
    maxiter::Int    # maximum number of iterations.
    verbose::Bool   # whether to show procedural information. 
    tol::Float64    # convergence tolerance
    theta::Float64  # switch parameter

    

end 
    
