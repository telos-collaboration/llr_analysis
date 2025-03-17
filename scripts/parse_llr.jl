using Pkg; Pkg.activate(".")
using LLRParsing

path = "/media/fabian/Adata2TB/LLR/"
path = "/home/fabian/Documents/Physics/Data/DataLLR"

#h5fileSU3 = "output/SU3_llr_david.hdf5"
#h5fileSU3_sorted = "output/SU3_llr_david_sorted.hdf5"
#h5fileSp4 = "output/Sp4_llr_david.hdf5"
#h5fileSp4_sorted = "output/Sp4_llr_david_sorted.hdf5"
h5fileSp4_new = "output/Sp4_llr_new.hdf5"
h5fileSp4_sorted_new = "output/Sp4_llr_new_sorted.hdf5"

function llr_alldirs_hdf5(path,base_dir,file)
    for dir in readdir(joinpath(path,base_dir),join=true)
        llr_dir_hdf5(dir,file)
    end
end

skip_ens = nothing
#llr_alldirs_hdf5(path,"llr_su3",h5fileSU3)
#sort_by_central_energy_to_hdf5(h5fileSU3, h5fileSU3_sorted)
#llr_alldirs_hdf5(path,"llr_sp4_cleaned",h5fileSp4)
#sort_by_central_energy_to_hdf5(h5fileSp4, h5fileSp4_sorted)
#llr_alldirs_hdf5(path,"llr_sp4_new_cleaned",h5fileSp4_new)
isfile(h5fileSp4_sorted_new) && rm(h5fileSp4_sorted_new)
sort_by_central_energy_to_hdf5(h5fileSp4_new, h5fileSp4_sorted_new; skip_ens)


