n = rand(1:10)

A = rand(n, n)

r = polarfact(A, alg = :halley);


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

println("Halley method passed test...")

##########################################################

m = rand(1:10)

B = rand(m, m)

r2 = polarfact(B, alg = :qdwh);

# Test unitary matrix U

U2 = r2.U
H2 = r2.H

@test_approx_eq_eps U2'*U2 eye(m) 1e-7

# Test Hermitian positive semifefinite matrix H

@test issym(H2)
 
for i in eigvals(H2)
    @test i >= 0.
end

@test_approx_eq B U2*H2

println("QDWH method passed test...")
