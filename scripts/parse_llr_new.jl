using Pkg; Pkg.activate(".")
using LLRParsing

path   = "/media/fabian/Adata2TB/LLR_New/CSD3/"
h5file = "output/Sp4_llr_new.hdf5"
h5file_sorted = "output/Sp4_llr_new_sorted.hdf5"

for dir in readdir(path,join=true)
    llr_dir_hdf5(dir,h5file)
end
skip_ens=[
          "4x80_11repeats_228replicas",
          "5x64_20repeats_95replicas",
          "5x72_22repeats_152replicas",
          "5x80_20repeats_76replicas",
          "6x72_26repeats_48replicas",
          "6x72_20repeats_76replicas",
          "6x84_20repeats_76replicas",
          ]
sort_by_central_energy_to_hdf5(h5file,h5file_sorted;skip_ens)