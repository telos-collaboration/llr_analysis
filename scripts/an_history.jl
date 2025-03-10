using HDF5
using Plots
using Statistics
function a_vs_central_action(h5dset,run;ind=nothing)
    N_replicas = read(h5dset[run],"N_replicas")
    N_repeats = read(h5dset[run],"N_repeats")
    a_last = zeros(N_replicas,N_repeats)
    S_last = zeros(N_replicas,N_repeats)
    # read all last elements for a and the central action
    for i in 1:N_replicas, j in 1:N_repeats
        n_steps = length(h5dset[run]["$(j-1)/Rep_$(i-1)/is_rm"])
        ind = isnothing(ind) ? n_steps : min(ind,n_steps) 
        a_last[i,j] = h5dset[run]["$(j-1)/Rep_$(i-1)/a_sorted"][ind]
        S_last[i,j] = h5dset[run]["$(j-1)/Rep_$(i-1)/S0_sorted"][ind]
    end
    # average over all repeats
    S0  = S_last[:,1]
    a0  = dropdims(mean(a_last,dims=2),dims=2)
    Δa0 = dropdims(std(a_last,dims=2),dims=2)/sqrt(N_repeats)
    return a0, Δa0, S0, ind
end
function a_vs_central_action_plot(h5dset,runs;indices)
    plt = plot()
    for run in runs
        for ind in indices
            a0, Δa0, S0, ind = a_vs_central_action(h5dset,run;ind)
            Nt = read(fid[run],"Nt")
            Nl = read(fid[run],"Nl")
            V  = Nl^3 * Nt
            up = S0/(6V)
            plot!(plt,up,a0,yerr=Δa0,marker=:auto,label="$(Nt)x$(Nl): ΔE=$(round(2(S0[2]-S0[1])/6V,sigdigits=1)) (Nr+RM steps=$ind)")
        end
    end
    return plt
end
function a_vs_central_action_plot!(plt,h5dset,run;index=nothing)
    a0, Δa0, S0, ind = a_vs_central_action(h5dset,run;ind=index)
    Nt = read(h5dset[run],"Nt")
    Nl = read(h5dset[run],"Nl")
    V  = Nl^3 * Nt
    up = S0/(6V)
    plot!(plt,up,a0,yerr=Δa0,marker=:auto,label="$(Nt)x$(Nl): ΔE=$(round(2(S0[2]-S0[1])/6V,sigdigits=1)) (Nr+RM steps=$ind)")
end