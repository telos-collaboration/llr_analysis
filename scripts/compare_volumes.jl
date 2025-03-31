using LLRParsing
using HDF5
using Plots
gr(fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

file = "output/test_sorted.hdf5"
h5id = h5open(file)
runs = keys(h5id)

filter!(contains("5x"),runs)
#runs = runs[[2,4,5,6,8]] # Nt4
runs = runs[[1,5,6,7]]   # Nt5

plt = plot()
for run in runs
    use_lens = false
    LLRParsing.a_vs_central_action_plot!(plt,h5id,run;index=nothing,highlight_index=nothing,lens=use_lens)
end

plot!(plt,xlims=(0.569,0.574) ,ylims=(7.333,7.345),legend=:bottomright) # Nt4
plot!(plt,xlims=(0.5888,0.590),ylims=(7.487,7.492),legend=:bottomright) # Nt5