#
# Newton's method
#
# Reference:
# Nicholas J. Higham, Computing the Polar Decomposition ---with Applications,
# SIAM J. Sci. Statist. Comput. Vol. 7, Num 4 (1986) pp. 1160-1174.
#
   
function newtion!{T}(X::Matrix{T}, # size (n,n)
                     U::Matrix{T}, # size (n,n)
                     H::Matrix{T}, # size (n,n)
                     maxiter::Int, # the maximum number of iterations
                     iter::Int,    # the number of iterations
                     tol::{T},  
                     verbose::Bool
                     )
