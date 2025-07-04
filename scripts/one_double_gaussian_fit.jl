using LLRParsing
using HDF5
using Plots
using ArgParse
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function plot_all_histogram_fits(file, plotfile, run; fit, A1=1, A2=1, name="fit")
    plotdir = dirname(plotfile)
    ispath(plotdir) || mkpath(plotdir)
    fid  = h5open(file)
    try
        plt = plot(title=LLRParsing.fancy_title(run))
        βc  = LLRParsing.beta_at_equal_heights(fid,run;A1,A2)
        plot_plaquette_histogram!(plt,fid,run,βc)
        if fit 
            ups, f, Δf = LLRParsing.histogram_jackknife_fit(fid,run,βc)
            LLRParsing.plot_double_gaussian_fit_difference(plt,fid,run, βc)
            plot!(plt,ups,f,ribbon=Δf,label="double Gaussian fit",lw=2)
        end
        savefig(plotfile)
    catch
        @warn "Error for plot with peak-ratio $A1:$A2 for $run"
    end
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
        "--run_name"
            help = "Run to be analysed"
            required = true
        "--relative_peak_height"
            help = "Ratio of peak heights used for finding the critical coupling"
            arg_type = Int
            default = 1
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    plotfile = args["plot_file"]
    run_name = args["run_name"]
    A1 = args["relative_peak_height"]
    plot_all_histogram_fits(file, plotfile, run_name; fit=true, A1=A1, A2=1, name="")
end
main()