using LLRParsing
using HDF5
using Plots
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function plot_all_histogram_fits(fid,runs)
    for run in runs
        plt = plot()
        βc  = LLRParsing.beta_at_equal_heights(fid,run;A1=1,A2=1)
        plot_plaquette_histogram!(plt,fid,run,βc;xlims=(0.5885,0.5905))
        LLRParsing.plot_double_gaussian_fit!(plt,fid,run,βc)
        display(plt)
    end

    for run in runs[2:end]
        plt = plot()
        βc  = LLRParsing.beta_at_equal_heights(fid,run;A1=2,A2=1)
        plot_plaquette_histogram!(plt,fid,run,βc;xlims=(0.5885,0.5905))
        LLRParsing.plot_double_gaussian_fit!(plt,fid,run,βc)
        display(plt)
    end
end

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
runs  = keys(fid)
plot_all_histogram_fits(fid,runs)