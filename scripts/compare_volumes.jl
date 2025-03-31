using LLRParsing
using HDF5
using Plots
gr(fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

file = "output/test_sorted.hdf5"
h5id = h5open(file)
runs = keys(h5id)
filter!(contains("4x"),runs)

plt = plot()
for run in runs
    use_lens = run == last(runs)
    LLRParsing.a_vs_central_action_plot!(plt,h5id,run;index=nothing,highlight_index=nothing,lens=use_lens)
end
plt