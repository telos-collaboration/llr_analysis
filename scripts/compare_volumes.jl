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

file = "output/test_sorted.hdf5"
h5id = h5open(file)
runs = keys(h5id)

runsNt4 = filter(contains("4x"),runs)
runsNt5 = filter(contains("5x"),runs)
runsNt6 = filter(contains("6x"),runs)

pltNt4 = a_vs_central_action_plot(h5id,runsNt4,lens=false)
pltNt5 = a_vs_central_action_plot(h5id,runsNt5,lens=false)
pltNt6 = a_vs_central_action_plot(h5id,runsNt6,lens=false)

plot!(pltNt4,xlims=(0.569,0.574) ,ylims=(7.333,7.345),legend=:bottomright) # Nt4
plot!(pltNt5,xlims=(0.5888,0.590),ylims=(7.487,7.492),legend=:bottomright) # Nt5
plot!(pltNt6,legend=:bottomright) # pltNt6

display(pltNt4)
display(pltNt5)
display(pltNt6)
