using Pkg; Pkg.activate(".")
using HDF5
using LLRParsing
using BenchmarkTools
using ProgressMeter
using HiRepParsing

function _all_files_from_dict(dir,replica_dirs)
    files = AbstractString[]
    for repeat in keys(replica_dirs), rep in replica_dirs[repeat]
        push!(files, joinpath(dir,repeat,rep,"out_0"))
    end
    return files
end
function parse_entire_llr_dir_to_hdf5(dir,h5file;suffix="") 
    fid = h5open(h5file,"cw")

    # get all repeats and replicas and store that information for future use
    repeats, replica_dirs = get_repeat_and_replica_dirs(dir)
    files      = _all_files_from_dict(dir,replica_dirs)
    N_repeats  = length(repeats)
    # assure the global lattice parameters are identical for all repeats and replicas
    N_replicas = only(unique([length(replica_dirs[r]) for r in repeats]))
    Nt = only(unique(first.(latticesize.(files))))  
    Nl = only(unique(last.(latticesize.(files))))

    name = "$(Nt)x$(Nl)_$(N_repeats)repeats_$(N_replicas)replicas"*suffix
    write(fid,joinpath(name,"N_repeats"),N_repeats)
    write(fid,joinpath(name,"N_replicas"),N_replicas)
    write(fid,joinpath(name,"Nt"),Nt)
    write(fid,joinpath(name,"Nl"),Nl)

    @showprogress desc="parsing $name" for repeat in repeats
        for rep in replica_dirs[repeat]
            file = joinpath(dir, repeat,rep,"out_0")
            dS0, S0, plaq, a, is_rm, S0_fxa, a_fxa, poly = parse_llr(file)
            write(fid,joinpath(name,repeat,rep,"dS0"),dS0)
            write(fid,joinpath(name,repeat,rep,"S0"),S0)
            write(fid,joinpath(name,repeat,rep,"plaq"),plaq)
            write(fid,joinpath(name,repeat,rep,"a"),a)
            write(fid,joinpath(name,repeat,rep,"is_rm"),is_rm)
            write(fid,joinpath(name,repeat,rep,"S0_fxa"),S0_fxa)
            write(fid,joinpath(name,repeat,rep,"a_fxa"),a_fxa)
            write(fid,joinpath(name,repeat,rep,"poly"),poly)
        end
    end
    close(fid)
end
function test_parsing(dir)
    repeats, replica_dirs = get_repeat_and_replica_dirs(dir)
    @showprogress for repeat in repeats
        for rep in replica_dirs[repeat]
            file = joinpath(dir, repeat,rep,"out_0")
            ans1 = parse_llr(file)
            ans2 = parse_llr_full(file)
            @assert ans1 == ans2
        end
    end
    
end

path = "/home/fabian/Documents/Physics/Data/DataLLR"
h5fileSp4 = "output/Sp4_llr_david.hdf5"
h5fileSU3 = "output/SU3_llr_david.hdf5"
for dir in readdir(joinpath(path,"Sp4/sp4_backup_full"),join=true)
    parse_entire_llr_dir_to_hdf5(dir,h5fileSp4)
end
for dir in readdir(joinpath(path,"SU3/su3_backup_logfiles"),join=true)
    parse_entire_llr_dir_to_hdf5(dir,h5fileSU3)
end
