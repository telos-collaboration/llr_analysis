using Pkg; Pkg.activate(".")
using LLRParsing

path = "/home/fabian/Documents/Physics/Data/DataLLR"
path = "/media/fabian/Adata2TB/LLR/"

h5fileSp4 = "output/Sp4_llr_david.hdf5"
h5fileSU3 = "output/SU3_llr_david.hdf5"

for dir in readdir(joinpath(path,"llr_sp4"),join=true)
    llr_dir_hdf5(dir,h5fileSp4)
end
for dir in readdir(joinpath(path,"llr_su3"),join=true)
    llr_dir_hdf5(dir,h5fileSU3)
end
