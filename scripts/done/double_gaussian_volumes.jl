using LLRParsing
using HDF5
using Plots
using ArgParse
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function plot_all_histogram_fits(file, plotfile, β0, βmin, βmax; xmin, xmax)
    ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
    fid   = h5open(file)
    runs  = keys(fid)
    plt = plot(title="Comparison of probability distribution across volumes")
    for run in runs
        try
            βc  = LLRParsing.beta_at_equal_heights(fid,run,β0,βmin,βmax)
            plot_plaquette_histogram!(plt,fid,run,βc;xlims=(xmin, xmax))
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
        "--xmin"
            help = "lower limit for the plots x-axis"
            arg_type = Float64
            required = true
        "--xmax"
            help = "upper limit for the plots x-axis"
            arg_type = Float64
            required = true
        "--beta_min"
            help = "lower limit of β for finding the critical coupling"
            arg_type = Float64
            required = true
        "--beta_max"
            help = "upper limit of β for finding the critical coupling"
            arg_type = Float64
            required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    βmin, βmax = args["beta_min"], args["beta_max"]
    xmin, xmax = args["xmin"], args["xmax"]
    file    = args["h5file"]
    plotfile = args["plotfile"]
    β0 = (βmax+βmin)/2

    plot_all_histogram_fits(file, plotfile, β0, βmin, βmax; xmin, xmax)
end
main()