using Pkg; Pkg.activate(".")
using LLRParsing

path = "/media/fabian/Adata2TB/LLR_New/"
h5fileSp4 = "output/Sp4_llr_new.hdf5"

for dir in readdir(path,join=true)
    llr_dir_hdf5(dir,h5fileSp4)
end
