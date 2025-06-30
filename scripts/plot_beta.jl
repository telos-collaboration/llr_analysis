using DelimitedFiles
using Plots
using ArgParse
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=10,labelfontsize=12,left_margin=1Plots.mm)

function read_critical_betas(file)
    data  = readdlm(file,',',skipstart=1)
    Nt    = data[:,2]
    Nl    = data[:,3]
    ratio = data[:,4]
    βc    = data[:,6]
    Δβc   = data[:,7]
    return only(unique(Nt)), Nl, βc, Δβc, only(unique(ratio))
end
function plot_critical_beta!(plt,file)
    @show file
    Nt, L, βc, Δβc, ratio = read_critical_betas(file)
    tks = (inv.(L), (L"1/%$Li" for Li in L))
    scatter!(plt,inv.(L),βc,xticks=tks,yerr=Δβc,markershape=:hexagon,label=L"N_t = %$Nt:$P_\beta$ peak ratio %$ratio:1")
end
function plot_critical_beta(files,plotfile)
    ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
    plt = plot(ylabel=L"critical $\beta_c$",xlabel=L"inverse spatial volume $1/N_l$")
    for file in files
        plot_critical_beta!(plt,file)
    end
    plot!(plt;xflip=true)
    savefig(plt,plotfile)
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--plotfile"
            help = "Where to save the plot"
            required = true
        "arg"
            help = "CSV file(s) containing the values of critical beta"
            required = true
            nargs = '+'
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    plotfile = args["plotfile"]
    file     = args["arg"]
    plot_critical_beta(file,plotfile)
end
main()