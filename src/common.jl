# common types or functions

# the result type

immutable Result{T}
    U::Matrix{T}
    H::Matrix{T}
    niters::Int
    converged::Bool
end
