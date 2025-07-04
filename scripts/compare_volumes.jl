using LLRParsing
using LaTeXStrings
using HDF5
using Plots
using ArgParse
gr(size=(425,282),fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=10,legendfontsize=7,tickfontsize=7,labelfontsize=10,left_margin=0Plots.mm)

a_vs_central_action_plot(h5id,runs::Vector;kws...) = a_vs_central_action_plot!(plot(),h5id,runs;kws...)
function a_vs_central_action_plot!(plt,h5id,runs::Vector;kws...)
    # default plot limits
    xmin, xmax = +Inf, -Inf 
    ymin, ymax = +Inf, -Inf 
    for run in runs
        a0, Δa0, S0, _ = a_vs_central_action(h5id,run)
        Nt = read(h5id[run],"Nt")
        Nl = read(h5id[run],"Nl")
        Nrep = read(h5id[run],"N_replicas")
        up = @. S0/(6*Nl^3*Nt)
        label=L"N_l\!=\!%$Nl,N_{\!\mathrm{rep}}\!=\!%$Nrep"
        LLRParsing.a_vs_central_action_plot!(plt,a0, Δa0, S0, Nt, Nl, Nrep; label, kws...)
        # find useful plot limits for the volume comparison
        ind  = findmax(Δa0)[2]
        δind = min(ind,length(a0)-ind)÷4
        xmin, xmax = min(xmin,up[ind-δind]), max(xmax,up[ind+δind])
        ymin, ymax = min(ymin,a0[ind-δind]), max(ymax,a0[ind+δind])
    end
    plot!(plt,xlims=(xmin,xmax),ylims=(ymin,ymax))
    return plt
end     
function an_action_volumes(file,plotdest;title)
    ispath(dirname(plotdest)) || mkpath(dirname(plotdest))
    h5id = h5open(file)
    runs = keys(h5id)
    plt  = a_vs_central_action_plot(h5id,runs,lens=false)
    title = latexstring(title)
    plot!(plt;legend=:bottomright,xlabel=L"u_p",ylabel=L"a_n",title)
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
    end
    return parse_args(s)
end
function main()
    args    = parse_commandline()
    file    = args["h5file"]
    plotdst = args["plot_file"]
    title   = args["title"]
    an_action_volumes(file,plotdst;title)
end
main()