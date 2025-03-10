using Pkg; Pkg.activate(".")
using LLRParsing

path   = "/media/fabian/Adata2TB/LLR_New/CSD3/"
h5file = "output/Sp4_llr_new.hdf5"
h5file_sorted = "output/Sp4_llr_new_sorted.hdf5"

for dir in readdir(path,join=true)
    llr_dir_hdf5(dir,h5file)
end
sort_by_central_energy_to_hdf5(h5file,h5file_sorted)
