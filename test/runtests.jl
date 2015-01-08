using PolarFact
using Base.Test

include("test_newton.jl")
include("test_halley.jl")
include("test_svd.jl")

println("Success!")
