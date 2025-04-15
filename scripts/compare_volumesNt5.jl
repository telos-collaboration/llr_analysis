using LLRParsing
using HDF5
using Plots
gr(fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

a_vs_central_action_plot(h5id,runs::Vector;kws...) = a_vs_central_action_plot!(plot(),h5id,runs;kws...)
function a_vs_central_action_plot!(plt,h5id,runs::Vector;kws...)
    for run in runs
        LLRParsing.a_vs_central_action_plot!(plt,h5id,run;kws...)
    end
    return plt
end     

isdir("plots") || mkdir("plots")
file = "output/test_sorted.hdf5"
h5id = h5open(file)
runs = keys(h5id)

runsNt5 = filter(contains("5x"),runs)
pltNt5 = a_vs_central_action_plot(h5id,runsNt5,lens=false)
plot!(pltNt5,xlims=(0.5887,0.59005),ylims=(7.487,7.492),legend=:bottomright) # Nt5
savefig(pltNt5,joinpath("plots","Nt5.pdf"))
display(pltNt5)