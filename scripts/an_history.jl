using HDF5
using Plots
using Statistics
gr(frame=:box)
function a_vs_central_action(h5dset,run)
    N_replicas = read(h5dset[run],"N_replicas")
    N_repeats = read(h5dset[run],"N_repeats")
    a_last = zeros(N_replicas,N_repeats)
    S_last = zeros(N_replicas,N_repeats)
    # read all last elements for a and the central action
    for i in 1:N_replicas, j in 1:N_repeats
        a_last[i,j] = h5dset[run]["$(j-1)/Rep_$(i-1)/a"][end]
        S_last[i,j] = h5dset[run]["$(j-1)/Rep_$(i-1)/S0"][end]
    end
    # Sort by the central action to account for different swaps
    # over different repeats
    for j in 1:N_repeats
        perm = sortperm(S_last[:,j])
        S_last[:,j] = S_last[perm,j]
        a_last[:,j] = a_last[perm,j]
    end
    # make sure that the sorted central action alwas matches
    for i in 1:N_replicas
        @assert allequal(S_last[i,:])
    end
    # average over all repeats
    S0  = S_last[:,1]
    a0  = dropdims(mean(a_last,dims=2),dims=2)
    Δa0 = dropdims(std(a_last,dims=2),dims=2)/sqrt(N_repeats)
    return a0, Δa0, S0
end

file = "output/SU3_llr_david.hdf5"
fid  = h5open(file)
runs = keys(fid)

plt = plot()
for run in runs[1:1]
    a0, Δa0, S0 = a_vs_central_action(fid,run)
    Nt = read(fid[run],"Nt")
    Nl = read(fid[run],"Nl")
    V  = Nl^3 * Nt
    up = S0/(6V)
    plot!(plt,up,a0,yerr=Δa0,marker=:auto,label="$(Nt)x$(Nl): ΔE=$(round(2(S0[2]-S0[1])/6V,sigdigits=1))")
end
plt