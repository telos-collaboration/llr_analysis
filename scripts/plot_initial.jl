using LLRParsing
using Plots
using HDF5
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

h5file_out = "output/test_sorted.hdf5"
h5dset = h5open(h5file_out)
runs = keys(h5dset)

plt = LLRParsing.a_vs_central_action_plot!(plot(),h5dset,runs[4];lens=false,index=1,highlight_index=nothing)
display(plt)

function read_initial_a(h5dset,run;repeat)
    replicas = keys(h5dset[run][repeat])
    Nl = h5dset[runs[1]]["Nl"][]
    Nt = h5dset[runs[1]]["Nt"][]
    a0 = [ read(h5dset[run][repeat][repl],"a0") for repl in replicas]
    S0 = [ h5dset[run][repeat][repl]["S0_sorted"][1] for repl in replicas]
    u0 = S0/(6*Nl^3*Nt)
    return a0, u0
end

repeat = "0"
a0, u0 = read_initial_a(h5dset,runs[4];repeat)
scatter(u0,a0)
u0