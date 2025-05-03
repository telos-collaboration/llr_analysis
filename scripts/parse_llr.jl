using Pkg; Pkg.activate(".")
using LLRParsing
using HiRepOutputCleaner

path = "/home/fabian/Documents/Physics/Data/DataCSD/CSD3/"
path = "/home/fabian/Documents/Physics/Data/DataLLR/llr_sp4_cleaned"
path = "/home/fabian/Downloads/LLR_cont"
#newpath = "./output/LLRout"
#clean_llr_directory(path,newpath;checkpoint_pattern=nothing,last_ranges=nothing,warn=false)
#path    = newpath

isdir("./output") || mkdir("./output")
h5file        = "output/test_Nt5.hdf5"
h5file_sorted = "output/test_Nt5_sorted.hdf5"

function llr_alldirs_hdf5(path,base_dir,file)
    for dir in readdir(joinpath(path,base_dir),join=true)
        if occursin("5x",basename(dir))
            @show basename(dir)
            llr_dir_hdf5(dir,file)
        end
    end
end

isfile(h5file) && rm(h5file)
isfile(h5file_sorted) && rm(h5file_sorted)
llr_alldirs_hdf5(path,"",h5file)
sort_by_central_energy_to_hdf5(h5file, h5file_sorted)
#include("trajectory_overview.jl")
#include("compare_volumesNt5.jl")