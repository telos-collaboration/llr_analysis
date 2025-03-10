using HDF5
using Plots
using LaTeXStrings
using Statistics
gr(fontfamily="Computer Modern",frame=:box,titlefontsize=11)

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
    plot!(plt,a[repeat+1],label="",lw=1,markersize=1,marker=:circle,msw=0.1,title=run,xlabel="NR/RM iteration m",ylabel=L"a_n^{(m)}")
    return plt
end
function plot_a_trajectory_all!(plt,h5dset,run,replica)
    a = a_trajectory(h5dset,run;replica)
    plot!(plt, a,label="",lw=1,markersize=1,marker=:circle,msw=0.1,title=run,xlabel="NR/RM iteration m",ylabel=L"a_n^{(m)}")
    return plt
end
function plot_nr_rm_shading!(plt,h5dset,run,repeat,replica)
    isrm = read(h5dset[run],"$repeat/Rep_$replica/is_rm")
    nr, rm = findfirst(isrm), length(isrm)
    vspan!(plt,[1,nr+1], color = :green, alpha = 0.2, labels = "NR");
    vspan!(plt,[nr+1,rm],color = :blue,  alpha = 0.2, labels = "RM");
end
function plot_a_repeat_average!(plt,h5dset,run;replica)
    a         = a_trajectory(h5dset,run;replica)
    a_last    = last.(a)
    N_repeats = length(a_last) 
    a0, Δa    = mean(a_last), std(a_last)/sqrt(N_repeats)
    hspan!(plt,[a0-Δa,a0+Δa], color = :gray, alpha = 0.4, labels = L"a_n");
end

h5file_out = "output/Sp4_llr_david_sorted.hdf5"
h5dset = h5open(h5file_out)
runs = keys(h5dset)
run  = runs[end]

# Plot 
replica_id = 0
repeat_id  = 4

plt = plot()
plot_a_trajectory_repeat!(plt,h5dset,run,repeat_id,replica_id)
plot_nr_rm_shading!(plt,h5dset,run,repeat_id,replica_id)
plot_a_repeat_average!(plt,h5dset,run;replica=replica_id)