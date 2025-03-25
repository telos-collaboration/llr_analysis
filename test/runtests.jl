using DelimitedFiles
using LLRParsing
using Test
using CodecZstd

@testset "All tests" begin
    include("std_tests.jl")
    include("llr_tests.jl")
    include("test_rm_filtering.jl")
end