using Pkg; Pkg.activate(".")
using LLRParsing

path = "/home/fabian/Documents/Physics/Data/DataLLR"
path = "/media/fabian/Adata2TB/LLR/"

h5file  = "data_assets/Sp4_std_data.hdf5"
basedir = joinpath(path,"std_sp4")
importance_sampling_dir_hdf5(basedir,h5file)

h5file  = "data_assets/SU3_std_data.hdf5"
basedir = joinpath(path,"std_su3")
importance_sampling_dir_hdf5(basedir,h5file)