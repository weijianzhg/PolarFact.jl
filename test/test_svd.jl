n = rand(1:10)
m = n + rand(0:5)
A = rand(m, n)

r = polarfact(A, alg = :svd);


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

println("SVD method passed test...")
