using Pkg; Pkg.activate(".")
using HDF5
using LLRParsing

h5file = "LLR_data.hdf5"
path = "/home/fabian/Documents/Physics/"
basedir = joinpath(path,"Data/DataLLR/ImportanceSampling/Importance_Sampling_logs/")
importance_sampling_dir_hdf5(basedir,h5file)