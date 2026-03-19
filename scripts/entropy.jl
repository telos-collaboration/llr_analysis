using ArgParse
using Plots
using LLRParsing
using DelimitedFiles
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
    palette = :Set1_5,
)
LINESTYLES = [:solid, :dot, :dash]
MARKERS = [:circle, :diamond, :dtriangle, :heptagon, :hexagon, :ltriangle, :octagon, :pentagon, :rect, :rtriangle, :star4, :star5, :star6, :star7, :star8, :utriangle ]

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--plot_file"
        help = "Where to save the plots"
        required = true
        "--critical_entropy"
        help = "CSV containing entropy at the critical point"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    plotfile = args["plot_file"]
    entropy_csv = args["critical_entropy"]
    critical_entropy = Float64(readdlm(entropy_csv, skipstart = 1)[1, 1])
    return LLRParsing.plot_entropy(file, plotfile, critical_entropy)
end
main()
