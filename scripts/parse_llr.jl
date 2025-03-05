using Test
using NaturalSort
using Parsers
function get_repeat_and_replica_dirs(base_dir)
    # Obtain all directories containing repeats and sort them
    repeat_dirs = filter( str -> all(isdigit, str), readdir(base_dir))
    sort!(repeat_dirs,lt=natural)
    
    dir_dict = Dict{String,Vector{String}}()
    for repeat in repeat_dirs
        rx = r"Rep_[0-9]+"
        repeat_path  = joinpath(base_dir,repeat)  
        replica_dirs = filter(startswith(rx),readdir(repeat_path))
        sort!(replica_dirs,lt=natural) 
        dir_dict[repeat] = replica_dirs 
    end
    return repeat_dirs, dir_dict
end
function parse_dS0(file)
    dS0 = NaN
    pattern = "[MAIN][0]LLR Delta S"
    for line in eachline(file)
        if startswith(line,pattern)
            dS0 = parse(Float64,line[length(pattern)+1:end])
            return dS0
        end
    end
    return dS0
end
function parse_S0(file)
    pattern = "[SWAP][10]New Rep Par S0 = "
    pos1    = length(pattern)
    S0      = Float64[]
    is_fxa = false
    for line in eachline(file)
        if startswith(line,"[SYSTEM][0]Process finalized.")
            is_fxa = false
        end
        if startswith(line,"[MAIN][0]Robins Monro update done.")
            is_fxa = true
        end        
        if !is_fxa && startswith(line,pattern)
            pos2 = first(findnext("dS",line,pos1))
            append!(S0,parse(Float64,line[pos1:pos2-1]))
        end
    end
    return S0
end
function parse_NR_plaquette(file)
    pattern = "[MAIN][0]NR Plaq a fixed "
    pos     = length(pattern)
    plaq    = Float64[]
    for line in eachline(file)
        if startswith(line,pattern)
            append!(plaq,parse(Float64,line[pos:end]))
        end
    end
    return plaq
end
function parse_RM_plaquette(file)
    pattern = "[MAIN][0]Plaq a fixed "
    pos     = length(pattern)
    plaq    = Float64[]
    for line in eachline(file)
        if startswith(line,pattern)
            append!(plaq,parse(Float64,line[pos:end]))
        end
    end
    return plaq
end
function parse_a_NR(file)
    pattern = "[MAIN][0]NR <a_rho"
    a       = Float64[]
    for line in eachline(file)
        if startswith(line,pattern)
            pos = findfirst('=',line) + 1
            append!(a,parse(Float64,line[pos:end]))
        end
    end
    return a
end
function parse_a_RM(file)
    pattern = "[MAIN][0]<a_rho"
    a       = Float64[]
    for line in eachline(file)
        if startswith(line,pattern)
            pos = findfirst('=',line) + 1
            append!(a,parse(Float64,line[pos:end]))
        end
    end
    return a
end
function parse_fixeda_S0_a_dS(file;S0_last,a_last,dS_last)
    is_fxa = false
    S0 = Float64[S0_last]
    a  = Float64[a_last]
    dS = Float64[dS_last]
    rx = r" S0 ([0-9]+.[0-9]+),  a  ([0-9]+.[0-9]+) , dS ([0-9]+.[0-9]+)"
    for line in eachline(file)
        if startswith(line,"[SYSTEM][0]Process finalized.")
            is_fxa = false
        end
        if startswith(line,"[MAIN][0]Robins Monro update done.")
            is_fxa = true
        end
        if is_fxa && startswith(line,"[llr:setreplica][0]New LLR Param:")
            vals = match(rx,line).captures 
            append!(S0,parse(Float64,vals[1]))
            append!(a ,parse(Float64,vals[2]))
            append!(dS,parse(Float64,vals[3]))
         end
    end
    # remove the last entry, because no measurement has been performed after the last fixed-a update
    return S0[1:end-1], a[1:end-1], dS[1:end-1]
end
function _parse_data!(array,string;n)
    opts = Parsers.Options(delim=' ', ignorerepeated=true)
    io = IOBuffer(string)
    for i in 1:n
        array[i] = Parsers.parse(Float64, io, opts)
    end
end
function parse_fun_polyakov_loop(file)
    poly    = ComplexF64[]
    pattern = "[FUND_POLYAKOV][0]Polyakov direction 0 = "
    pos     = length(pattern)
    vals    = zeros(2)
    for line in eachline(file)
        if startswith(line,pattern)
            _parse_data!(vals,line[pos:end];n=2)
            append!(poly, vals[1] + im*vals[2])
        end
    end
    return poly
end

base_dir = "/home/fabian/Downloads/llr_parser_test_data/su3_4x20_8/"
repeats, replica_dirs = get_repeat_and_replica_dirs(base_dir)

# raw logs
fileSU3 = "/home/fabian/Downloads/llr_parser_test_data/su3_4x20_8/0/Rep_0/out_0"
fileSp4 = "/home/fabian/Downloads/llr_parser_test_data/sp4_4x20_48/0/Rep_0/out_0"

# David's comparison files
david_fa_Sp4 = "test/test_data/4x20_0_Rep0_fa.csv"
david_RM_Sp4 = "test/test_data/4x20_0_Rep0_RM.csv"

using DelimitedFiles
fa = readdlm(david_fa_Sp4,',',skipstart=1)
rm = readdlm(david_RM_Sp4,',',skipstart=1)
# only use rows where the replica (the 5th or 7th entry respectively) number is zero
fa = fa[fa[:,5] .== 0, :]
rm = rm[rm[:,7] .== 0, :]
rm = sortslices(rm,dims=1,by=x->(x[1]))
fa = sortslices(fa,dims=1,by=x->(x[8]))

dS0 = parse_dS0(fileSp4)
S0  = parse_S0(fileSp4)
plaq_RM = parse_RM_plaquette(fileSp4)
plaq_NR = parse_NR_plaquette(fileSp4)
plaq = vcat(plaq_NR,plaq_RM)
a_NR = parse_a_NR(fileSp4)
a_RM = parse_a_RM(fileSp4)
a    = vcat(a_NR,a_RM) 
poly = parse_fun_polyakov_loop(fileSp4)
S0_fxa, a_fxa, dS_fxa = parse_fixeda_S0_a_dS(fileSp4;S0_last=S0[end],a_last=a[end],dS_last=dS0)
a_fxa_david = Float64[]
S0_fxa_david = Float64[]
for i in  Iterators.partition(eachindex(fa[:,1]),100) 
    append!(a_fxa_david,  only(unique(fa[i,1])))
    append!(S0_fxa_david,  only(unique(fa[i,2])))
end

@testset "Compare LLR parsing with David's code" begin  
    @test dS0 == 2only(unique(rm[:,4]))
    @test S0 == rm[:,3]
    @test isapprox(plaq,rm[:,6]./rm[:,5]./6)
    @test a == -rm[:,2]
    @test abs.(poly) == fa[:,6]
    @test only(unique(dS_fxa)) == 2only(unique(rm[:,4]))
    @test a_fxa == -a_fxa_david
    @test S0_fxa == S0_fxa_david
end