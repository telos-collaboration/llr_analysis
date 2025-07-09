using ArgParse
using Plots
using LLRParsing
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 1Plots.mm,
)

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--plot_file"
        help = "Where to save the plot"
        required = true
        "--run_name"
        help = "Run to be analysed"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    plotfile = args["plot_file"]
    run = args["run_name"]
    return LLRParsing.plot_free_energy(file, plotfile, run)
end
main()
