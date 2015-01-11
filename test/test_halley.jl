n = rand(1:10)

A = rand(n, n)

r = polarfact(A, alg = :halley);


# Test unitary matrix U

U = r.U

@test_approx_eq_eps U'*U eye(n) 1e-7

# Test Hermitian positive semifefinite matrix H

@test issym(r.H)
 
for i in eigvals(r.H)
    @test i >= 0.
end

println("Halley method passed test...")
