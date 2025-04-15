using HDF5
using LLRParsing

file = "output/SU3_llr_sorted.hdf5"
fid  = h5open(file)

run     = "4x20_20repeats_55replicas"
repeat  = "0"
replica = "Rep_0"

# I have checked that those are the correct values that David's code is using
a_all, S_all = LLRParsing.a_vs_central_action_repeats(fid,run;ind=nothing)[1:2]
a = a_all[:,end]
S = S_all[:,end]
