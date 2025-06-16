using Pkg; Pkg.activate(".")
using LLRParsing
using HDF5
using Plots
using ArgParse
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function plot_all_histogram_fits(file, plotdir, β0, βmin, βmax; xmin, xmax, fit, A1=1, A2=1, name="fit")
    ispath(plotdir) || mkpath(plotdir)
    fid   = h5open(file)
    runs  = keys(fid)
    for run in runs
        try
            plt = plot(title=LLRParsing.fancy_title(run))
            βc  = LLRParsing.beta_at_equal_heights(fid,run,β0,βmin,βmax;A1,A2)
            plot_plaquette_histogram!(plt,fid,run,βc;xlims=(xmin, xmax))
            if fit 
                ups, f, Δf = LLRParsing.histogram_jackknife_fit(fid,run,βc)
                LLRParsing.plot_double_gaussian_fit_difference(plt,fid,run, βc)
                plot!(plt,ups,f,ribbon=Δf,label="double Gaussian fit",lw=2)
            end
            savefig(joinpath(plotdir,"$(name)_$run.pdf"))
        catch
            @warn "Error for plot $name for $run"
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
            help = "Where to save the plot"
            required = true
        "--xmin"
            help = "lower limit for the plots x-axis"
            arg_type = Float64
            required = true
        "--xmax"
            help = "upper limit for the plots x-axis"
            arg_type = Float64
            required = true
        "--beta_min"
            help = "lower limit of β for finding the critical coupling"
            arg_type = Float64
            required = true
        "--beta_max"
            help = "upper limit of β for finding the critical coupling"
            arg_type = Float64
            required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    βmin, βmax = args["beta_min"], args["beta_max"]
    xmin, xmax = args["xmin"], args["xmax"]
    file    = args["h5file"]
    plotdir = args["plot_dir"]
    β0 = (βmax+βmin)/2

    plot_all_histogram_fits(file, plotdir, β0, βmin, βmax; xmin, xmax, fit=true, name = "fit")
    plot_all_histogram_fits(file, plotdir, β0, βmin, βmax; xmin, xmax, fit=false, A1=2, A2=1, name="two-to-one")
end
main()