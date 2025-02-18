using Pkg; Pkg.activate(".")
using BenchmarkTools
using DelimitedFiles
using LLRParsing
using Test

compare_file = "/home/fabian/Documents/Physics/Analysis/LLR_SU3/LLRAnalysis/data/IS/IS_4x20/7.32/full.csv"
file = "/home/fabian/Documents/Physics/Data/DataLLR/ImportanceSampling/Importance_Sampling_noCSV_noRW/IS_4x20/7.32/output_file"

Nt, Nl, plaq, beta, poly = parse_importance_sampling(file)
poly1 = abs.(real.(poly))

data = readdlm(compare_file,',',skipstart=1)
plaq2 = data[:,2]
poly2 = data[:,3]

Δpoly = maximum(abs.(poly1 - poly2))
Δplaq = maximum(abs.(plaq - plaq2))

@testset "Compare with David's files" begin
    @test iszero(Δplaq)
    @test maximum(Δpoly) < 1E-8
end