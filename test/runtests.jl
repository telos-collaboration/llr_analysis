using DelimitedFiles
using LLRParsing
using Test
using CodecZstd

# Test parsing of importance sampling files 
compare_file = "./test_data/importance_sampling_4x20_beta7.32_full.csv"
file         = "./test_data/importance_sampling_4x20_beta7.32_logs.zst"
file_tmp     = "./test_data/tmp.txt"
compare_file_std = "./test_data/4x20_std.csv"

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

# Now test the "standard observables" from David's file ''
ans = std_observables(plaq, abs.(real.(poly)), Nt, Nl) 
binder_plaq, Δbinder_plaq, sh_plaq, Δsh_plaq, plaq_vev, Δplaq_vev, poly_sus, Δpoly_sus, poly_vev, Δpoly_vev = ans
data_std = readdlm(compare_file_std,',',skipstart=1)

@testset "Compare results with David's" begin
    @testset "Raw Polyakov loop and plaquette" begin
        @test iszero(Δplaq)
        @test maximum(Δpoly) < 1E-8
    end
    @testset "'std_observables:'" begin
        @test isapprox(plaq_vev,data_std[1,2]) 
        @test isapprox(Δplaq_vev,data_std[1,3],rtol=0.3) 
        @test isapprox(sh_plaq,data_std[1,4],rtol=1E-5) 
        @test isapprox(Δsh_plaq,data_std[1,5],rtol=0.5) 
        @test isapprox(poly_vev,data_std[1,6]) 
        @test isapprox(Δpoly_vev,data_std[1,7],rtol=5E-2) 
        @test isapprox(poly_sus,data_std[1,8],rtol=1E-5)
        @test isapprox(Δpoly_sus,data_std[1,9],rtol=1E-2) 
        @test isapprox(binder_plaq,data_std[1,10])
        @test isapprox(Δbinder_plaq,data_std[1,11],rtol=0.5)
    end
end