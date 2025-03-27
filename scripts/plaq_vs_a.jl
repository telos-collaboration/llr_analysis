using HDF5
using LLRParsing
using DelimitedFiles

h5file_out = "output/test_sorted.hdf5"
h5dset = h5open(h5file_out,"r")

for run in keys(h5dset)
    @show run
    a0, Δa0, S0, ind = a_vs_central_action(h5dset,run)
    Nt = read(h5dset[run],"Nt")
    Nl = read(h5dset[run],"Nl")
    V  = Nl^3 * Nt
    up = S0/(6V)

    @assert issorted(up)
    data = hcat(a0,up)
    writedlm("output/$run.txt",data,',')
end

