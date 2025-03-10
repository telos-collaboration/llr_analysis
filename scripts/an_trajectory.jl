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

using HDF5
using Plots


h5file_out = "output/Sp4_llr_david_sorted.hdf5"
h5dset = h5open(h5file_out)
runs = keys(h5dset)

# Plot 
replica = 0
run  = runs[1]
a    = a_trajectory(h5dset,run;replica)
isrm = read(h5dset[run],"0/Rep_$replica/is_rm")
plt  = plot(a,label="",lw=1,markersize=1,marker=:circle,msw=0.1,title=run,xlabel="NR/RM iteration m",ylabel="an^(m)")
nr, rm = findfirst(isrm), length(isrm)
vspan!(plt,[1,nr+1], color = :green, alpha = 0.2, labels = "NR");
vspan!(plt,[nr+1,rm],color = :blue,  alpha = 0.2, labels = "RM");
plt