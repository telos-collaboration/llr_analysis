using HDF5
using Plots
using LaTeXStrings
using Statistics
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)
include("an_history.jl")
function a_trajectory(h5dset,run;replica=0)
    N_replicas = read(h5dset[run],"N_replicas")
    N_repeats = read(h5dset[run],"N_repeats")
    @assert replica < N_replicas
    a = [ Float64[] for i in 1:N_repeats]
    for j in 1:N_repeats
        a[j] = vec(h5dset[run]["$(j-1)/Rep_$(replica)/a_sorted"][])
    end
    return a
end
function plot_a_trajectory_repeat!(plt,h5dset,run,repeat,replica)
    a    = a_trajectory(h5dset,run;replica)
    plot!(plt,a[repeat+1],label="",lw=1,markersize=1,marker=:circle,msw=0.1)
    return plt
end
function plot_a_trajectory_all!(plt,h5dset,run,replica)
    a = a_trajectory(h5dset,run;replica)
    plot!(plt, a,label="",lw=1,markersize=1,marker=:circle,msw=0.1)
    return plt
end
function plot_nr_rm_shading!(plt,h5dset,run,repeat,replica)
    isrm = read(h5dset[run],"$repeat/Rep_$replica/is_rm")
    nr, rm = findfirst(isrm), length(isrm)
    vspan!(plt,[1,nr+1], color = :green, alpha = 0.2, labels = "NR");
    vspan!(plt,[nr+1,rm],color = :blue,  alpha = 0.2, labels = "RM");
    return plt
end
function plot_a_repeat_average!(plt,h5dset,run;replica)
    a         = a_trajectory(h5dset,run;replica)
    a_last    = last.(a)
    N_repeats = length(a_last) 
    a0, Δa    = mean(a_last), std(a_last)/sqrt(N_repeats)
    hspan!(plt,[a0-Δa,a0+Δa], color = :black, alpha = 0.8, labels = L"a_n");
    return plt
end
function fancy_title(run)
    rx = r"([0-9]+)x([0-9]+)_([0-9]+)repeats_([0-9]+)replicas"
    m  = match(rx,run).captures
    Lt, Ls, Nrep, Nint = m
    return L"%$Lt \times %$(Ls)^3: N_{\mathrm{intervals}}=%$Nint, N_{\mathrm{rep}}=%$Nrep"
end
function full_trajectory_plot(h5dset,run,repeat_id,replica_id;lens=true)
    plt1 = plot_a_trajectory_repeat!(plot(),h5dset,run,repeat_id,replica_id)
    plt2 = plot_a_trajectory_all!(plot(),h5dset,run,replica_id)
    plt3 = a_vs_central_action_plot!(plot(),h5dset,run;lens,index=nothing,highlight_index=replica_id)
    plt4 = a_variance_vs_central_action_plot!(plot(),h5dset,run;index=nothing,highlight_index=replica_id)
    for plt in [plt1,plt2]
        plot_nr_rm_shading!(plt,h5dset,run,repeat_id,replica_id)
        plot_a_repeat_average!(plt,h5dset,run;replica=replica_id)
    end
    plot!(plt1, ylims=ylims(plt2), title=fancy_title(run),ylabel=L"a_n^{(m)}")
    plot!(plt2, ylims=ylims(plt2), title="", xlabel="NR/RM iteration m",ylabel=L"a_n^{(m)}")  
    plot!(plt3, title="", xlabel="",ylabel=L"a_n")  
    plot!(plt4, title="", xlabel=L"\textrm{central}~\textrm{plaquette}~u_p")
    if lens
        plot!(plt3, subplot=2, left_margin=0Plots.mm, tickfontsize=7, ylabel="")  
    end
    l   =  grid(4, 1, heights=[0.29 ,0.29, 0.29, 0.13])
    plt = plot(plt1, plt2, plt3, plt4, layout = l, size=(500,1000))
    return plt
end

h5file_out = "output/SU3_llr_david_sorted.hdf5"
h5file_out = "output/Sp4_llr_david_sorted.hdf5"
h5file_out = "output/Sp4_llr_new_sorted.hdf5"

h5dset = h5open(h5file_out)
runs = keys(h5dset)
run  = runs[1]

# Plot 
Nreplicas        = read(h5dset[run],"N_replicas")
a0, Δa0, S0, ind = a_vs_central_action(h5dset,run)
replica_id       = min(findmax(Δa0)[2],Nreplicas-1)
repeat_id        = 0
plt1 = full_trajectory_plot(h5dset,run,repeat_id,replica_id,lens=true)
plt2 = full_trajectory_plot(h5dset,run,repeat_id,Nreplicas-1)
plt3 = full_trajectory_plot(h5dset,run,repeat_id,1)
plot!(plt2,legend=:topleft)
plt = plot(plt3, plt1, plt2, layout = grid(1, 3), size=(1400,1000))