using DelimitedFiles
using Plots
using ArgParse
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=10,labelfontsize=12,left_margin=1Plots.mm)

function read_critical_betas(file)
    data  = readdlm(file,',',skipstart=1)
    N     = size(data,1)
    Nl    = data[1:N÷2,2]
    βc11  = data[1:N÷2,5]
    Δβc11 = data[1:N÷2,6]
    βc22  = data[N÷2+1:end,5]
    Δβc22 = data[N÷2+1:end,6]
    return Nl, βc11, Δβc11, βc22, Δβc22
end
function plot_critical_beta(file,plotfile,Nt)
    L, βc11, Δβc11, βc22, Δβc22 = read_critical_betas(file)
    plt = plot(;title=L"N_t = %$Nt",ylabel=L"critical $\beta_c$",xlabel=L"inverse spatial volume $1/N_l$")
    tks = (inv.(L), (L"1/%$Li" for Li in L))
    scatter!(plt,inv.(L),βc11,xticks=tks,yerr=Δβc11,markershape=:hexagon,label=L"$P_\beta$ peak ratio 1:1")
    scatter!(plt,inv.(L),βc22,xticks=tks,yerr=Δβc22,markershape=:circle,label=L"$P_\beta$ peak ratio 2:1")
    plot!(plt;xflip=true)
    ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
    savefig(plt,plotfile)
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--file"
            help = "CSV file containing the values of critical beta"
            required = true
        "--plotfile"
            help = "Where to save the plot"
            required = true
        "--Nt"
            help = "Temporal lattice extent"
            arg_type = Int
            required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    Nt       = args["Nt"]
    file     = args["file"]
    plotfile = args["plotfile"]
    plot_critical_beta(file,plotfile,Nt)
end
main()