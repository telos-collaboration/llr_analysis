using DelimitedFiles
using LLRParsing
using Test
using CodecZstd

# Test parsing of importance sampling files 
compare_file = "./test_data/importance_sampling_4x20_beta7.32_full.csv"
file         = "./test_data/importance_sampling_4x20_beta7.32_logs.zst"
include("std_tests.jl")