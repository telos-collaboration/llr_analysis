using HDF5
using Plots
using LaTeXStrings
using Statistics
gr(fontfamily="Computer Modern",frame=:box,titlefontsize=12,legendfontsize=11,labelfontsize=13)

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
    hspan!(plt,[a0-Δa,a0+Δa], color = :gray, alpha = 0.8, labels = L"a_n");
    return plt
end
function fancy_title(run)
    rx = r"([0-9]+)x([0-9]+)_([0-9]+)repeats_([0-9]+)replicas"
    m  = match(rx,run).captures
    Lt, Ls, Nrep, Nint = m
    return L"%$Lt \times %$(Ls)^3: N_{\mathrm{int}}=%$Nint, N_{\mathrm{rep}}=%$Nrep"
end
function full_trajectory_plot(h5dset,run,repeat_id,replica_id)
    plt1 = plot_a_trajectory_repeat!(plot(),h5dset,run,repeat_id,replica_id)
    plt2 = plot_a_trajectory_all!(plot(),h5dset,run,replica_id)
    for plt in [plt1,plt2]
        plot_nr_rm_shading!(plt,h5dset,run,repeat_id,replica_id)
        plot_a_repeat_average!(plt,h5dset,run;replica=replica_id)
    end
    plot!(plt1, ylims=ylims(plt2), title=fancy_title(run),ylabel=L"a_n^{(m)}")
    plot!(plt2, ylims=ylims(plt2), title="", xlabel="NR/RM iteration m",ylabel=L"a_n^{(m)}")  
    l   = @layout [a ; b]
    plt = plot(plt1, plt2, layout = l, legend = :outerright)
    return plt
end

h5file_out = "output/Sp4_llr_david_sorted.hdf5"
h5dset = h5open(h5file_out)
runs = keys(h5dset)
run  = runs[end]

# Plot 
Nreplicas  = read(h5dset[run],"N_replicas")
repeat_id  = 0
replica_id = Nreplicas ÷2
plt = full_trajectory_plot(h5dset,run,repeat_id,replica_id)
