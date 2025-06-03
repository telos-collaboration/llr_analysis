using Pkg; Pkg.activate(".")
using HDF5
using LLRParsing
using Plots
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function histogram_comparison(fid,runs,betas)
    plt = plot()
    for (i,run) in enumerate(runs)
        xl = (0.5885,0.5905)
        plot_plaquette_histogram!(plt,fid,run,betas[i];xlims=xl)
    end
    return plt
end

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
betas = [7.48967, 7.48970, 7.48982, 7.48969, 7.48975]
runs  = keys(fid)
plt = histogram_comparison(fid,runs,betas)
savefig("Nt5_histogram.pdf")
display(plt)