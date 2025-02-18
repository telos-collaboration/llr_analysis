using Pkg; Pkg.activate(".")
using HDF5
using LLRParsing

function parse_importance_sampling_hdf5(basedir,h5file;outname="output_file")
    isfile(h5file) && rm(h5file)
    # Iterate recursively through all directories in basedir
    # Check if there is a matching output file
    # If yes: Save data to the hdf5 file 
    for (root,dirs,files) in walkdir(basedir)
        if outname ∈ files
            infile = joinpath(root,outname)
            importance_sampling_hdf5(infile,h5file)
        end
    end
end

using BenchmarkTools
basedir = "/home/fabian/Documents/Physics/Data/DataLLR/ImportanceSampling/Importance_Sampling_noCSV_noRW/"
h5file = "test.hdf5"
@profview parse_importance_sampling_hdf5(basedir,h5file)
@time parse_importance_sampling_hdf5(basedir,h5file)
f = h5open(h5file)
close(f)