using LLRParsing
using HDF5
using Plots
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function main()
    file = "data_assets/test_Nt5_sorted.hdf5"
    fid  = h5open(file)
    runs = keys(fid)
    plt  = plot()
    xl   = (0.5885,0.5905)
    for run in runs
        βc  = LLRParsing.beta_at_equal_heights(fid,run)
        plot_plaquette_histogram!(plt,fid,run,βc, xlims=xl)
    end
    display(plt)
end
main()
