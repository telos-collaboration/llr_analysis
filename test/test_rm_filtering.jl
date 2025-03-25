# Use two files that is known to be corrupted
# file1 can be sufficiently cleaned by HiRepOutputCleaner.jl 
# file2 has some residual corruption in the Robbins Monro iteration

# These tests clean both files and trhen assure that the number of correctly measured steps matches
# even though the initial corrupted files do not.


file1 = "test_data/LLR_5x72_152_0_Rep0_out_0.txt"
file2 = "test_data/LLR_5x72_152_0_Rep2_out_0.txt"
S0_Rep0 = parse_llr_corrupted(file1)[2]
S0_Rep2 = parse_llr_corrupted(file2)[2]

@testset begin
    @test length(S0_Rep0) == 58
    @test length(S0_Rep2) == 58
end
