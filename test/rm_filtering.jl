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
    @testset "Parse with data corruption in Robbins Monro update" begin
        @test length(S0_Rep0) == 58
        @test length(S0_Rep2) == 58
    end
    @testset "Check for regression on known uncorrupted file" begin 
        file = "./test_data/4x20_0_Rep0_out_0.txt"
        ans1 = parse_llr(file)
        ans2 = parse_llr_corrupted(file)
        @test ans1 == ans2
    end
    @testset "Check identification of incomplete LLR update in logfile" begin
        file = "./test_data/LLR_4x80_228_0_Rep1_out_0.txt"
        invalid_lines = LLRParsing.find_invalid_llr_update_lines(file)
        @test invalid_lines == [825, 826, 827]
    end
end