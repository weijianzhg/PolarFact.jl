n = rand(1:10)

A = rand(n,n)

r = polarfact(A);

# Test unitary matrix U

U = r.U

@test_approx_eq_eps U'*U eye(n) 1e-3

# Test Hermitian positive semifefinite matrix H

@test issym(r.H)
 
for i in eigvals(r.H)
    @test i >= 0.
end

println("Newton method passed test...")
