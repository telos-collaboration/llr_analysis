using LLRParsing
using Plots
using HDF5
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=7Plots.mm)

h5file_out = "data_assets/test_Nt5_sorted.hdf5"
h5dset = h5open(h5file_out)
runs = [last(keys(h5dset))]

function overview(h5dset,run;repeat_id=1, replica_id = read(h5dset[run],"N_replicas")÷2)
    Nreplicas  = read(h5dset[run],"N_replicas")
    plt1 = full_trajectory_plot(h5dset,run,repeat_id,replica_id ,lens=false)
    plt2 = full_trajectory_plot(h5dset,run,repeat_id,Nreplicas-1,lens=false)
    plt3 = full_trajectory_plot(h5dset,run,repeat_id,1          ,lens=false)
    plot!(plt2,legend=:topleft)
    plt = plot(plt3, plt1, plt2, layout = grid(1, 3), size=(1400,1000))
    return plt
end
for run in runs
   plt = overview(h5dset,run,repeat_id=4,replica_id=24)
   display(plt)
   savefig("$(run)_overview.svg")
end