using DelimitedFiles
using LLRParsing
using Test
using CodecZstd
using HDF5

@testset "All tests" begin
    include("std_tests.jl")
    include("llr_tests.jl")
    include("rm_filtering.jl")
    include("replica_filtering.jl")
end