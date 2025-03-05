using Pkg; Pkg.activate(".")
using HDF5
using LLRParsing

h5file  = "output/Sp4_std_data.hdf5"
basedir = "/home/fabian/Dokumente/Physics/Data/DataLLR/ImportanceSampling/Importance_Sampling_noCSV_noRW/"
importance_sampling_dir_hdf5(basedir,h5file)

h5file  = "output/SU3_std_data.hdf5"
basedir = "/media/fabian/Adata2TB/std"
importance_sampling_dir_hdf5(basedir,h5file)