using ArgParse
using Plots
using LLRParsing
gr(fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=1Plots.mm)

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
    file     = args["h5file"]
    plotfile = args["plot_file"]
    run      = args["run_name"]
    LLRParsing.plot_free_energy(file, plotfile, run)
end
main()