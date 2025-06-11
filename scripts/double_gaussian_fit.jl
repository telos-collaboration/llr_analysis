using LLRParsing
using HDF5
using Plots
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function plot_all_histogram_fits(file, β0, βmin, βmax; xl)
    fid   = h5open(file)
    runs  = keys(fid)
    for run in runs
        @show run
        plt = plot()
        βc  = LLRParsing.beta_at_equal_heights(fid,run,β0,βmin,βmax;A1=1,A2=1)
        ups, f, Δf = LLRParsing.histogram_jackknife_fit(fid,run,βc)
        plot_plaquette_histogram!(plt,fid,run,βc;xlims=xl)
        LLRParsing.plot_double_gaussian_fit_difference(plt,fid,run, βc)
        plot!(plt,ups,f,ribbon=Δf,label="double Gaussian fit",lw=2)
        display(plt)
    end
end

β0, βmin, βmax = 7.490, 7.488, 7.492
file  = "data_assets/Sp4_Nt5_sorted.hdf5"
plot_all_histogram_fits(file, β0, βmin, βmax; xl=(0.5885,0.5905))

#β0, βmin, βmax = 7.337, 7.340, 7.343
#file  = "data_assets/Sp4_Nt4_sorted.hdf5"
#plot_all_histogram_fits(file, β0, βmin, βmax, xl=(0.568,0.576))