using Pkg; Pkg.activate(".")
using LLRParsing

path = "/media/fabian/Adata2TB/LLR/"
path = "/home/fabian/Documents/Physics/Data/DataLLR/"

h5fileSU3        = "output/SU3_llr.hdf5"
h5fileSU3_sorted = "output/SU3_llr_sorted.hdf5"

function llr_alldirs_hdf5(path,base_dir,file)
    for dir in readdir(joinpath(path,base_dir),join=true)
        llr_dir_hdf5(dir,file)
    end
end

llr_alldirs_hdf5(path,"llr_su3",h5fileSU3)
sort_by_central_energy_to_hdf5(h5fileSU3, h5fileSU3_sorted)