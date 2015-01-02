# PolarFact

A Julia package for polar decomposition.


## Types

* Result::

	immutable Result
		U::Matrix{Float64}
		H::Matrix{Float64}
		niters::Int
		converged::Bool
	end

## Examples

```julia
julia> using MatrixDepot

julia> using PolarFact

julia> A = matrixdepot("moler", 7)
7x7 Array{Float64,2}:
  1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 -1.0   2.0   0.0   0.0   0.0   0.0   0.0
 -1.0   0.0   3.0   1.0   1.0   1.0   1.0
 -1.0   0.0   1.0   4.0   2.0   2.0   2.0
 -1.0   0.0   1.0   2.0   5.0   3.0   3.0
 -1.0   0.0   1.0   2.0   3.0   6.0   4.0
 -1.0   0.0   1.0   2.0   3.0   4.0   7.0

julia> r = polarfact(A);

julia> r.U
7x7 Array{Float64,2}:
 1.0          1.82835e-16  9.15431e-17  …  2.35503e-17  1.28443e-17  8.56205e-18
 1.82835e-16  1.0          4.57841e-17     1.17784e-17  6.42389e-18  4.2822e-18 
 9.15431e-17  4.57841e-17  1.0             5.89728e-18  3.21635e-18  2.14404e-18
 4.60348e-17  2.30237e-17  1.15277e-17     2.9656e-18   1.61743e-18  1.07819e-18
 2.35503e-17  1.17784e-17  5.89728e-18     1.0          8.27437e-19  5.51574e-19
 1.28443e-17  6.42389e-18  3.21635e-18  …  8.27437e-19  1.0          3.00826e-19
 8.56205e-18  4.2822e-18   2.14404e-18     5.51574e-19  3.00826e-19  1.0  

julia> r.niters
15
```



