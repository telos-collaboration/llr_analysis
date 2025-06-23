using Pkg; Pkg.activate(".")
using HDF5
using LLRParsing
using Plots
using LaTeXStrings
using PDFmerger
using ProgressMeter
using ArgParse
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function plot_repeat(run,repeat,Nrep;kws...)
    title = LLRParsing.fancy_title(run)*", repeat #$repeat"
    a0    = hcat([read(fid[run],"$repeat/Rep_$i/a_sorted") for i in 0:Nrep-1]...)
    NRpRM = first(size(a0))
    plt   = plot(a0;ms=1,legend=:outerright,title,label="",xlabel="updates (NR + RM = $NRpRM)",ylabel=L"a_n")
    plot!(plt;kws...)
end
function plot_all_an_trajectories(file,outdir;kws...)
    fid  = h5open(file)
    @showprogress "Plot a_n sequences" for run in keys(fid)
        outfile = joinpath(outdir,"$run.pdf")
        ispath(outdir) || mkpath(outdir)
        isfile(outfile) && rm(outfile)
        
        repeats = read(fid[run],"repeats")
        Nrep = read(fid[run],"N_replicas")
        
        for r in repeats
            plt = plot_repeat(run,r,Nrep;kws...)
            savefig(plt, "tmp.pdf")
            append_pdf!(outfile, "tmp.pdf", cleanup=true)
        end
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
    plot_all_an_trajectories(args["h5file"],args["plot_dir"])
end
main()