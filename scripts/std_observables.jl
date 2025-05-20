using Pkg; Pkg.activate(".")
using LLRParsing

h5file     = "data_assets/SU3_std_data.hdf5"
h5file_out = "data_assets/SU3_std_obs.hdf5"
LLRParsing.write_all_std_observables_to_hdf5(h5file,h5file_out)

h5file     = "data_assets/Sp4_std_data.hdf5"
h5file_out = "data_assets/Sp4_std_obs.hdf5"
LLRParsing.write_all_std_observables_to_hdf5(h5file,h5file_out)