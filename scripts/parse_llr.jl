using Pkg; Pkg.activate(".")
using LLRParsing
using HiRepOutputCleaner

path    = "/home/fabian/Documents/Physics/Data/DataCSD/CSD3/"
newpath = "./output/LLRout" 
clean_llr_directory(path,newpath;checkpoint_pattern=nothing,last_ranges=nothing,warn=false)
path    = newpath

h5fileSp4_new        = "output/test.hdf5"
h5fileSp4_sorted_new = "output/test_sorted.hdf5"

function llr_alldirs_hdf5(path,base_dir,file)
    for dir in readdir(joinpath(path,base_dir),join=true)
        llr_dir_hdf5(dir,file)
    end
end

isfile(h5fileSp4_new) && rm(h5fileSp4_new)
llr_alldirs_hdf5(path,"",h5fileSp4_new)
isfile(h5fileSp4_sorted_new) && rm(h5fileSp4_sorted_new)
sort_by_central_energy_to_hdf5(h5fileSp4_new, h5fileSp4_sorted_new)
include("trajectory_overview.jl")