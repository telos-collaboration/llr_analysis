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

base_dir = "/home/fabian/Downloads/llr_parser_test_data/su3_4x20_8/"
repeats, replica_dirs = get_repeat_and_replica_dirs(base_dir)

fileSU3 = "/home/fabian/Downloads/llr_parser_test_data/su3_4x20_8/0/Rep_0/out_0"
fileSp4 = "/home/fabian/Downloads/llr_parser_test_data/sp4_4x20_48/0/Rep_0/out_0"
@test parse_dS0(fileSU3) == 1216.86282
@test parse_dS0(fileSp4) == 122.55319

