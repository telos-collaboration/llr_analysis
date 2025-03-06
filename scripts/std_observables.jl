using Pkg; Pkg.activate(".")
using LLRParsing

h5file     = "output/SU3_std_data.hdf5"
h5file_out = "output/SU3_std_obs.hdf5"
LLRParsing.write_all_std_observables_to_hdf5(h5file,h5file_out)

h5file     = "output/Sp4_std_data.hdf5"
h5file_out = "output/Sp4_std_obs.hdf5"
LLRParsing.write_all_std_observables_to_hdf5(h5file,h5file_out)