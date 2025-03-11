using Pkg; Pkg.activate(".")
using LLRParsing

path = "/media/fabian/Adata2TB/LLR/"
path = "/home/fabian/Documents/Physics/Data/DataLLR"

h5fileSp4 = "output/Sp4_llr_david.hdf5"
h5fileSU3 = "output/SU3_llr_david.hdf5"
h5fileSp4_new = "output/Sp4_llr_new.hdf5"
h5fileSp4_sorted = "output/SU3_llr_david_sorted.hdf5"
h5fileSp4_sorted = "output/SU3_llr_david_sorted.hdf5"
h5fileSU3_sorted_new = "output/Sp4_llr_new.hdf5"

for dir in readdir(joinpath(path,"llr_sp4_cleaned"),join=true)
    llr_dir_hdf5(dir,h5fileSp4)
end
#for dir in readdir(joinpath(path,"llr_sp4_new_cleaned"),join=true)
#    llr_dir_hdf5(dir,h5fileSp4_new)
#end
#for dir in readdir(joinpath(path,"llr_su3"),join=true)
#    llr_dir_hdf5(dir,h5fileSU3)
#end

sort_by_central_energy_to_hdf5(h5fileSp4, h5fileSp4_sorted)
sort_by_central_energy_to_hdf5(h5fileSU3, h5fileSU3_sorted)
sort_by_central_energy_to_hdf5(h5fileSp4_new, h5fileSp4_sorted_new)
