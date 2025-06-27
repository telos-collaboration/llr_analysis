using LLRParsing
using LaTeXStrings
using HDF5
using Plots
using ArgParse
gr(fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

a_vs_central_action_plot(h5id,runs::Vector;kws...) = a_vs_central_action_plot!(plot(),h5id,runs;kws...)
function a_vs_central_action_plot!(plt,h5id,runs::Vector;kws...)
    for run in runs
        LLRParsing.a_vs_central_action_plot!(plt,h5id,run;kws...)
    end
    return plt
end     
function an_action_volumes(file,plotdest;title,xmin,xmax,ymin,ymax)
    ispath(dirname(plotdest)) || mkpath(dirname(plotdest))
    h5id = h5open(file)
    runs = keys(h5id)
    plt  = a_vs_central_action_plot(h5id,runs,lens=false)
    title = latexstring(title)
    plot!(plt,xlims=(xmin,xmax),ylims=(ymin,ymax);title)
    plot!(plt,legend=:bottomright,xlabel=L"u_p",ylabel=L"a_n")
    savefig(plt,plotdest)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
            help = "HDF5 file containing the sorted results"
            required = true
        "--plot_file"
            help = "Where to save the plot"
            required = true
        "--title"
            help = "Title of the plot"
            required = true
        "--xmin"
            help = "lower limit for the plots x-axis"
            arg_type = Float64
            required = true
        "--xmax"
            help = "upper limit for the plots x-axis"
            arg_type = Float64
            required = true
        "--ymin"
            help = "lower limit for the plots y-axis"
            arg_type = Float64
            required = true
        "--ymax"
            help = "upper limit for the plots y-axis"
            arg_type = Float64
            required = true
    end
    return parse_args(s)
end
function main()
    args    = parse_commandline()
    file    = args["h5file"]
    plotdst = args["plot_file"]
    title   = args["title"]
    xmin     = args["xmin"]
    xmax     = args["xmax"]
    ymin     = args["ymin"]
    ymax     = args["ymax"]
    an_action_volumes(file,plotdst;title,xmin,xmax,ymin,ymax)
end
main()