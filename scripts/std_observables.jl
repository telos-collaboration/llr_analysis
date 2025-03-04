using Pkg; Pkg.activate(".")
using LLRParsing

h5file     = "LLR_data.hdf5"
h5file_out = "LLR_obs.hdf5"
write_all_std_observables_to_hdf5(h5file,h5file_out)