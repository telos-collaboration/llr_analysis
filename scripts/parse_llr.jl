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
"""
    Parse the new central action for the replica.
    In contrast to Davids parsing code I do not distinguish
    between NR, RM and fixed-a updates here but will later 
    split the output into corresping pices.
"""
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

base_dir = "/home/fabian/Downloads/llr_parser_test_data/su3_4x20_8/"
repeats, replica_dirs = get_repeat_and_replica_dirs(base_dir)

fileSU3 = "/home/fabian/Downloads/llr_parser_test_data/su3_4x20_8/0/Rep_0/out_0"
fileSp4 = "/home/fabian/Downloads/llr_parser_test_data/sp4_4x20_48/0/Rep_0/out_0"
@test parse_dS0(fileSU3) == 1216.86282
@test parse_dS0(fileSp4) == 122.55319
S0 = parse_S0(fileSp4)

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
parse_RM_plaquette(fileSp4)
parse_RM_plaquette(fileSU3)
parse_NR_plaquette(fileSp4)
parse_NR_plaquette(fileSU3)