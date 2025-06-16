using HDF5
using LLRParsing
using Plots
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

file = "data_assets/Sp4_Nt5_sorted.hdf5"
fid  = h5open(file)
runs = keys(fid)
run  = runs[3]

function plot_all_an_trajectories(fid,run;yl=(7.475,7.503))
    repeats = read(fid[run],"repeats")
    Nrep = read(fid[run],"N_replicas")
    for r in repeats
        title = LLRParsing.fancy_title(run)*", repeat #$r"
        a0    = hcat([read(fid[run],"$r/Rep_$i/a_sorted") for i in 0:Nrep-1]...)
        NRpRM = first(size(a0))
        plt   = plot(a0;ms=1,legend=:outerright,title,label="",xlabel="updates (NR + RM = $NRpRM)",ylabel=L"a_n")
        plot!(plt,ylims=yl)
        display(plt)
    end
end
plot_all_an_trajectories(fid,run)