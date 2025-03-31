using Pkg; Pkg.activate(".")
using LLRParsing
using HiRepOutputCleaner

path = "/home/fabian/Documents/Physics/Data/DataCSD/CSD3/"
path = "/home/fabian/Documents/Physics/Data/DataLLR/llr_sp4_cleaned"
newpath = "./output/LLRout"
clean_llr_directory(path,newpath;checkpoint_pattern=nothing,last_ranges=nothing,warn=false)
path    = newpath

isdir("./output") || mkdir("./output")
h5file        = "output/test.hdf5"
h5file_sorted = "output/test_sorted.hdf5"

function llr_alldirs_hdf5(path,base_dir,file)
    for dir in readdir(joinpath(path,base_dir),join=true)
        @show basename(dir)
        llr_dir_hdf5(dir,file)
    end
end

isfile(h5file) && rm(h5file)
isfile(h5file_sorted) && rm(h5file_sorted)
llr_alldirs_hdf5(path,"",h5file)
sort_by_central_energy_to_hdf5(h5file, h5file_sorted)
include("trajectory_overview.jl")