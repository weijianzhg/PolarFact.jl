n = rand(1:10)

A = Array(Float32, n, n)
copy!(A, rand(n,n))

r1 = polarfact(A, alg =:newton);

r2 = polarfact(A, alg =:qdwh);

r3 = polarfact(A, alg =:halley);

r4 = polarfact(A, alg =:svd);

r5 = polarfact(A, alg =:hybrid);

@test_approx_eq r1.U r2.U
@test_approx_eq r2.U r3.U
@test_approx_eq r3.U r4.U
@test_approx_eq r4.U r5.U

println("Type Float32 passed test...")


