# PolarFact

[![Build Status](https://travis-ci.org/weijianzhang/PolarFact.jl.svg?branch=master)](https://travis-ci.org/weijianzhang/PolarFact.jl)
| Julia 0.3 [![PolarFact](http://pkg.julialang.org/badges/PolarFact_release.svg)](http://pkg.julialang.org/?pkg=PolarFact&ver=release)
| Julia 0.4 [![PolarFact](http://pkg.julialang.org/badges/PolarFact_nightly.svg)](http://pkg.julialang.org/?pkg=PolarFact&ver=nightly)

A Julia package for the matrix polar decomposition.

## Install

To install the release version, type

```julia
julia> Pkg.add("PolarFact")
```

## Overview 

Every ``m-by-n`` matrix ``A`` has a polar decomposition ``A=UH``,
where the ``m-by-n`` matrix ``U`` has orthonormal columns if ``m>n``
or orthonormal rows if ``m<n`` and the ``n-by-n`` matrix ``H`` is
Hermitian positive semidefinite. For a square matrix ``A``, ``H`` is
unique. If in addition, ``A`` is nonsingular, then ``H`` is positive
definite and ``U`` is unique.

The polar decomposition is closely related to the singular value
decomposition (SVD). In particular, if ``A = P*S*Q'`` is a singular
value decomposition of A, then ``U = P*Q'`` and ``H = Q*S*Q'`` are the
corresponding polar factors. The orthonormal polar factor ``U`` is the
nearest orthonormal matrix to ``A`` in the Frobenius norm [1, Sec. 8.1]. 

[1] Nicholas J. Higham, Functions of Matrices: Theory and Computation,
SIAM, Philadelphia, PA, USA, 2008.

This package provides the following algorithms for computing matrix
polar decomposition:

* (Scaled) Newton's method

	Reference:
	[2] Nicholas J. Higham, Computing the Polar Decomposition ---with Applications,
	SIAM J. Sci. Statist. Comput. Vol. 7, Num 4 (1986) pp. 1160-1174.
	
* the Newton Schulz method 
  
    This method can only apply to matrix ``A`` such that ``norm(A) < sqrt(3)``.

	Reference:
	[3] Günther Schulz, Iterative Berechnung der reziproken Matrix, Z. Angew.
	Math. Mech.,13:57-59, (1933) pp. 114, 181.

* a hybrid Newton method

	Start with (scaled) Newton's method and switch to Newton-Schulz method
	when convergence is guaranteed.

	Reference:
	[4] Nicholas J. Higham and Robert S. Schreiber, Fast Polar
	Decomposition of an arbitrary matrix, SIAM, J. Sci. Statist. Comput.
	Vol. 11, No. 4 (1990) pp. 648-655

* Halley's method

	Reference:
	[5] Y. Nakatsukasa, Z. Bai and F. Gygi, Optimizing Halley's iteration 
	for computing the matrix polar decomposition, SIAM, J. Mat. Anal. 
	Vol. 31, Num 5 (2010) pp. 2700-2720. 

* the QR-based Dynamically weighted Halley (QDWH) method [5]  

* the SVD method

### Comments on Usage

The scaled Newton iteration is a well known and effective method for
computing the polar decomposition. It converges quadratically and is
backward stable under the assumption that the matrix inverses are
computed in a mixed backward forward stable way[6]. The QDWH is a
cubic-rate convergent method.  It is backward stable under the
assumption that column pivoting and either row pivoting or row sorting
are used in the QR factorization [6].  Without scaling, both type of
methods can be slow when the matrix is ill-conditioned.

On many modern computers, matrix multiplication can be performed
very efficiently. The Newton Schulz method requires two matrix
multiplication while the (scaled) Newton method requires one matrix
inversion. Thus the hybrid Newton is more efficient if matrix
multiplication is 2 times faster than the matrix inversion [4].

Comparing to the SVD approach, the iterative algorithms are much more
efficient when the matrix is nearly unitary (as in applications, for
example, where a time-dependent matrix drifts from orthogonality due
to rounding errors or other errors).

[6] Yuji Nakatsukasa and Nicholas J. Higham, Backward stability of
iterations for computing the polar decomposition, SIAM, J.
Matrix Anal. Appl. Vol. 33, No. 2, pp. 460-479. 


## Interface

The package provides a high-level function ``polarfact``:

```julia
	polarfact(A; alg, maxiter, tol, verbose)
```

The meaning of the arguments:

- ``A`` : the input matrix of type ``Matrix{Float64}``.

- ``alg``: a symbol that indicates the factorization algorithm (default = ``:newton``).

	This argument accepts the following values:

	- ``:newton``: scaled Newton's method
	- ``:qdwh``: the QR-based Dynamically weighted Halley (QDWH) method
	- ``:halley``: Halley's method
	- ``:schulz``: the Newton Schulz method
	- ``:hybrid``: a hybrid Newton method 
	- ``:svd``: the SVD method

- ``maxiter``: maximum number of iterations (default = ``100``).

- ``tol`` :  tolerance (default = ``1.0e-6``).

- ``verbose`` : whether to show procedural information (default = ``false``), where
               ``Iter`` is the number of iterations, ``Rel. err.`` is equal to
			   ``norm(preU - U, Inf) / norm(preU, Inf)`` and ``Obj.`` is equal to
			   ``norm(I - U'*U, Inf)``. 

*Note:* ``maxiter``, ``tol`` and ``verbose`` are not used for the
SVD method.

The output has type ``PolarFact.Result``, which is defined as 

```
	immutable Result
		U::Matrix{Float64}               # unitary factor
		H::Matrix{Float64}               # Hermitian positive semidefinite factor
		niters::Union(Int, Nothing)      # number of iterations or Nothing
		converged::Union(Bool, Nothing)  # whether the algorithm converges or Nothing
	end
```

*Note:* ``niters`` and ``converged`` are of type ``Nothing`` for the
SVD method. 

## Examples

```julia
	julia> using PolarFact

	julia> A = rand(10,10);

	julia> r = polarfact(A);

	julia> r.U
	10x10 Array{Float64,2}:
	-0.489361   0.0601629   0.346164   …  -0.0213289   0.596928    0.279623 
	0.302874  -0.129649    0.080602      -0.241901    0.0372137  -0.435765 
	0.111824  -0.236326   -0.302027       0.019033   -0.0162169   0.456752 
	-0.314102   0.452703    0.0132996      0.0460329  -0.481377    0.358254 
	-0.168465  -0.123994   -0.0880474      0.449138    0.0860416  -0.0702016
	0.153207  -0.053861    0.5924     …   0.660055   -0.080868   -0.128266 
	0.391474   0.361302    0.25525       -0.015308   -0.294591    0.202959 
	0.352644  -0.31368     0.424935      -0.294412    0.111984    0.489584 
	0.452478   0.206963   -0.411828       0.419035    0.388085    0.249209 
	0.153067   0.654781    0.0909015     -0.19659     0.382659   -0.174877 

	julia> r.niters
	6

	julia> r = polarfact(A, alg = :qdwh, verbose = true);
	Iter.    Rel. err.        Obj.         
	1     2.337305e-01     2.940272e+01
    2     5.454641e-01     5.108731e+00
    3     3.655509e-01     3.338175e-01
    4     4.833229e-02     1.140381e-03
    5     1.788780e-04     6.675496e-11
    6     1.047445e-11     2.144118e-15
```

## Acknowledgements

The design of the package is inspired by [NMF.jl](https://github.com/JuliaStats/NMF.jl).


