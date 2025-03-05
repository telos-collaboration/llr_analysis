using Pkg; Pkg.activate(".")
using HDF5
using LLRParsing
using BenchmarkTools
using ProgressMeter

function parse_entire_llr_dir_to_hdf5(dir,h5file;name=basename(dir)) 
    isfile(h5file) && rm(h5file)
    fid = h5open(h5file,"w")

    repeats, replica_dirs = get_repeat_and_replica_dirs(dir)
    @showprogress for repeat in repeats
        for rep in replica_dirs[repeat]
            file = joinpath(dir, repeat,rep,"out_0")
            dS0, S0, plaq, a, is_rm, S0_fxa, a_fxa, dS_fxa, poly = parse_llr_quick(file)
            write(fid,joinpath(name,repeat,rep,"dS0"),dS0)
            write(fid,joinpath(name,repeat,rep,"S0"),S0)
            write(fid,joinpath(name,repeat,rep,"plaq"),plaq)
            write(fid,joinpath(name,repeat,rep,"a"),a)
            write(fid,joinpath(name,repeat,rep,"is_rm"),is_rm)
            write(fid,joinpath(name,repeat,rep,"S0_fxa"),S0_fxa)
            write(fid,joinpath(name,repeat,rep,"a_fxa"),a_fxa)
            write(fid,joinpath(name,repeat,rep,"dS_fxa"),dS_fxa)
            write(fid,joinpath(name,repeat,rep,"poly"),poly)
        end
    end
    close(fid)
end

dir = "/home/fabian/Downloads/llr_parser_test_data/sp4_4x20_48"
h5file = "test.hdf5"
parse_entire_llr_dir_to_hdf5(dir,h5file)
