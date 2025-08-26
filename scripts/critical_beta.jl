using LLRParsing
using HDF5
using Plots
using ArgParse
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topright,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
)

function all_critical_beta(file, outfile; A1 = 1, A2 = 1)
    fid = h5open(file)
    runs = keys(fid)
    runs = filter(!startswith("provenance"), runs)
    header = !isfile(outfile)
    io = open(outfile, "a")
    print_provenance_csv(io)
    header && println(io, "group,run,Nt,Nl,A1,A2,βc,Δβc")
    for run in runs
        L = read(fid[run], "Nl")
        T = read(fid[run], "Nt")
        try
            βc, Δβc = LLRParsing.βc_jackknife(fid, run; A1, A2)
            println(io, "Sp4,$run,$T,$L,$A1,$A2,$βc,$Δβc")
        catch
            @warn "Cannot determine critical coupling for run $run with peak ratio $A1:$A2"
            println(io, "Sp4,$run,$T,$L,$A1,$A2,NaN,NaN")
        end
    end
    return close(io)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--outfile"
        help = "Where to save results"
        required = true
        "--peak1"
        help = "Ratio of peak heights used for finding the critical coupling"
        arg_type = Int
        default = 1
        "--peak2"
        help = "Ratio of peak heights used for finding the critical coupling"
        arg_type = Int
        default = 1
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    A1 = args["peak1"]
    A2 = args["peak2"]
    outfile = args["outfile"]
    isfile(outfile) && rm(outfile)
    return all_critical_beta(file, outfile; A1 = A1, A2 = A2)
end
main()
