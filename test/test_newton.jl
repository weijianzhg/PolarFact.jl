n = rand(1:10)

A = rand(n,n)

r = polarfact(A, alg = :newton);

# Test unitary matrix U

U = r.U
H = r.H

@test_approx_eq_eps U'*U eye(n) 1e-7

# Test Hermitian positive semifefinite matrix H

@test issym(H)
 
for i in eigvals(H)
    @test i >= 0.
end

@test_approx_eq A U*H

println("Newton method passed test...")
