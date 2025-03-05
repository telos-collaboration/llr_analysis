using Test
using NaturalSort
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
    for line in eachline(file)
        if startswith(line,pattern)
            pos2 = first(findnext("dS",line,pos1))
            append!(S0,parse(Float64,line[pos1:pos2-1]))
        end
    end
    return S0
end
function parse_RM_plaquette(file)
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
function parse_NR_plaquette(file)
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
function parse_fixeda_S0_a_dS(file)
    is_fxa = false
    S0 = Float64[]
    a  = Float64[]
    dS = Float64[]
    rx = r" S0 ([0-9]+.[0-9]+),  a  ([0-9]+.[0-9]+) , dS ([0-9]+.[0-9]+)"

    for line in eachline(file)
        if startswith(line,"[SYSTEM][0]Process finalized.")
            is_fxa = false
            @show is_fxa
        end
        if startswith(line,"[MAIN][0]Robins Monro update done.")
            is_fxa = true
            @show is_fxa
        end
        if is_fxa && startswith(line,"[llr:setreplica][0]New LLR Param:")
            vals = match(rx,line).captures 
            append!(S0,parse(Float64,vals[1]))
            append!(a ,parse(Float64,vals[2]))
            append!(dS,parse(Float64,vals[3]))
         end
    end
    return S0, a, dS
end

base_dir = "/home/fabian/Downloads/llr_parser_test_data/su3_4x20_8/"
repeats, replica_dirs = get_repeat_and_replica_dirs(base_dir)

fileSU3 = "/home/fabian/Downloads/llr_parser_test_data/su3_4x20_8/0/Rep_0/out_0"
fileSp4 = "/home/fabian/Downloads/llr_parser_test_data/sp4_4x20_48/0/Rep_0/out_0"
@test parse_dS0(fileSU3) == 1216.86282
@test parse_dS0(fileSp4) == 122.55319
S0 = parse_S0(fileSp4)
parse_RM_plaquette(fileSp4)
parse_RM_plaquette(fileSU3)
parse_NR_plaquette(fileSp4)
parse_NR_plaquette(fileSU3)
parse_a_NR(fileSp4)
parse_a_NR(fileSU3)
parse_a_RM(fileSp4)
parse_a_RM(fileSU3)
parse_fixeda_S0_a_dS(fileSp4)[1]
parse_fixeda_S0_a_dS(fileSU3)[1]