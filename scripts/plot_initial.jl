using LLRParsing
using Plots
using HDF5
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function read_initial_a(h5dset,run;repeat)
    replicas = keys(h5dset[run][repeat])
    Nl = h5dset[run]["Nl"][]
    Nt = h5dset[run]["Nt"][]
    a0 = [ read(h5dset[run][repeat][repl],"a0") for repl in replicas]
    S0 = [ h5dset[run][repeat][repl]["S0_sorted"][1] for repl in replicas]
    u0 = S0/(6*Nl^3*Nt)
    @show Nl, Nt
    return a0, u0
end

h5file_out = "data_assets/test_Nt5_sorted.hdf5"
h5dset = h5open(h5file_out)
runs   = keys(h5dset)
run    = runs[2]
repeat = "0"

plt = LLRParsing.a_vs_central_action_plot!(plot(),h5dset,run;lens=false,index=nothing,highlight_index=nothing)
#LLRParsing.a_vs_central_action_plot!(plt,h5dset,run;lens=false,index=nothing,highlight_index=nothing)
a0, u0 = read_initial_a(h5dset,run;repeat)
scatter!(u0,a0,label=L"inital $a_n$")
savefig("initial_an_Vega.pdf")
plt