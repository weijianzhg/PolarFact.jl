using PolarFact
using Base.Test

include("test_newton.jl")
include("test_halley.jl")
include("test_svd.jl")
include("test_hybrid.jl")
include("test_f32.jl") # test Float32 case

println("Success!")
