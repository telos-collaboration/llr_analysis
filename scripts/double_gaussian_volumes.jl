using LLRParsing
using HDF5
using Plots
using ArgParse
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function plot_all_histogram_fits(file, plotfile)
    ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
    fid   = h5open(file)
    runs  = keys(fid)
    plt = plot(title="Comparison of probability distribution across volumes")
    for run in runs
        try
            βc  = LLRParsing.beta_at_equal_heights(fid,run)
            plot_plaquette_histogram!(plt,fid,run,βc)
        catch
            @warn "Cannot determine critical β for $run"
        end
    end
    savefig(plt,plotfile)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
            help = "HDF5 file containing the sorted results"
            required = true
        "--plotfile"
            help = "Where to save the plot"
            required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    plotfile = args["plotfile"]
    plot_all_histogram_fits(file, plotfile)
end
main()