using HDF5
using Plots
using Statistics
using LaTeXStrings
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
function _lens_location(a0, Δa0, up)
    # identify approximate midpoint by the maximum of the variance of a_n
    # then estimate peak and dip position by looking for 
    #   1) the maximum before the midpoint
    #   2) the minimum after the midpoint
    mid,  m_ind = findmax(Δa0)
    peak, p_ind = findmax(a0[1:m_ind+1])
    dip,  d_ind = findmin(a0[m_ind:end])
    delta_peak_dip = max((m_ind-p_ind),d_ind)
    ylms = [dip,peak] .+ [dip-peak,peak-dip] ./ 3
    xinds = [p_ind - 3delta_peak_dip÷2, m_ind + 7delta_peak_dip÷2]
    xlms = up[xinds]
    return ylms, xlms, [dip, peak], up[[p_ind, m_ind + d_ind]]
end
function a_vs_central_action_plot!(plt,h5dset,run;index=nothing,highlight_index=nothing,lens=true)
    a0, Δa0, S0, ind = a_vs_central_action(h5dset,run;ind=index)
    Nt = read(h5dset[run],"Nt")
    Nl = read(h5dset[run],"Nl")
    V  = Nl^3 * Nt
    up = S0/(6V)
    plot!(plt,up,a0,yerr=Δa0,marker=:auto,label="ΔE=$(round(2(S0[2]-S0[1])/6V,sigdigits=1))")
    if !isnothing(highlight_index)
        mid = up[highlight_index]
        del = (up[highlight_index+1] - up[highlight_index])/2
        vspan!(plt,[mid-del,mid+del], color = :green, alpha = 0.7, labels = "replica");
    end
    if lens
        ylms, xlms, yticks, xticks = _lens_location(a0, Δa0, up)
        xticks = round.(xticks,sigdigits=4)
        yticks = round.(yticks,sigdigits=5)
        lens!(plt, xlms, ylms,inset = (1, bbox(0.60, 0.50, 0.38, 0.38 )))
        plot!(;subplot=2, yticks, xticks)
    end
end
function a_variance_vs_central_action_plot!(plt,h5dset,run;index=nothing,highlight_index=nothing)
    a0, Δa0, S0, ind = a_vs_central_action(h5dset,run;ind=index)
    Nt = read(h5dset[run],"Nt")
    Nl = read(h5dset[run],"Nl")
    V  = Nl^3 * Nt
    up = S0/(6V)
    plot!(plt,up,zero(a0),ribbon=Δa0,label=L"\Delta a_n")
    if !isnothing(highlight_index)
        mid = up[highlight_index]
        del = (up[highlight_index+1] - up[highlight_index])/2
        vspan!(plt,[mid-del,mid+del], color = :green, alpha = 0.8, labels = "replica");
    end
end