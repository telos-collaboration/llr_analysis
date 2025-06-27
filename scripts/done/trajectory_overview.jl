using LLRParsing
using Plots
using HDF5
using ArgParse
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=7Plots.mm)

function overview(h5dset,run;repeat_id=1, replica_id = read(h5dset[run],"N_replicas")÷2)
    Nreplicas  = read(h5dset[run],"N_replicas")
    plt1 = full_trajectory_plot(h5dset,run,repeat_id,replica_id ,lens=false)
    plt2 = full_trajectory_plot(h5dset,run,repeat_id,Nreplicas-1,lens=false)
    plt3 = full_trajectory_plot(h5dset,run,repeat_id,1          ,lens=false)
    plot!(plt2,legend=:topleft)
    plt = plot(plt3, plt1, plt2, layout = grid(1, 3), size=(1400,1000))
    return plt
end
function all_overview_plots(file,plotdir)
    h5dset = h5open(file)
    runs   = keys(h5dset)
    ispath(plotdir) || mkpath(plotdir)
    for run in runs
        Δa0 = a_vs_central_action(h5dset,run)[2]
        ind = findmax(Δa0)[2]
        plt = overview(h5dset,run,repeat_id=1,replica_id=ind)
        savefig(plt,joinpath(plotdir,"$run.pdf"))
    end
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
            help = "HDF5 file containing the sorted results"
            required = true
        "--plot_dir"
            help = "Where to save the plots"
            required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    all_overview_plots(args["h5file"],args["plot_dir"])
end
main()