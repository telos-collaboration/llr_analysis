using DelimitedFiles
using LLRParsing
using Test
using CodecZstd

# Test parsing of importance sampling files 
compare_file = "./test_data/importance_sampling_4x20_beta7.32_full.csv"
file         = "./test_data/importance_sampling_4x20_beta7.32_logs.zst"
file_tmp     = "./test_data/tmp.txt"

# Decompress file used in testing
input  = ZstdDecompressorStream(open(file, "r"))
output = open(file_tmp, "w")
write(output, input)
close(input)
close(output)

Nt, Nl, plaq, beta, poly = parse_importance_sampling(file_tmp)
data = readdlm(compare_file,',',skipstart=1)
rm(file_tmp)

poly1 = abs.(real.(poly))
plaq2 = data[:,2]
poly2 = data[:,3]

Δpoly = maximum(abs.(poly1 - poly2))
Δplaq = maximum(abs.(plaq - plaq2))

@testset "Compare with David's files" begin
    @test iszero(Δplaq)
    @test maximum(Δpoly) < 1E-8
end