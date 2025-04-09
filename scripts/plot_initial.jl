using LLRParsing
using Plots
using HDF5
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

h5file_out = "output/test_sorted.hdf5"
h5dset = h5open(h5file_out)
runs = keys(h5dset)

plt = LLRParsing.a_vs_central_action_plot!(plot(),h5dset,runs[4];lens=false,index=1,highlight_index=nothing)
display(plt)    