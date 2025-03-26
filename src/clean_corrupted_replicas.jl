# Utilities for 
# 1) Finding all indices for unique elements in an array v  
# 2) Finding the values of all items of v that occur more than once
# 3) Finding the inidces of all items of v that occur more than once allowing for slight deviations
unique_indices(v) = unique(i -> v[i], eachindex(v))
duplicate_vals(v) = v[setdiff(eachindex(v),unique_indices(v))]
approx_duplicate_indices(v) = findall(i-> any(isapprox(v[i]),duplicate_vals(v)), eachindex(v))

function read_non_matching_trajectory(fid;key)
    replicas     = keys(fid)
    traj_lengths = [length(fid["$rep/$key"]) for rep in replicas]
    nmin, nmax   = extrema(traj_lengths)
    
    T   = eltype(fid["$(first(replicas))/$key"])
    res = zeros(T,(length(replicas),nmax)) .+ typemax(T)
    
    for (i,rep) in enumerate(replicas)
        val = read(fid,"$rep/$key")
        res[i,1:length(val)] = val
    end
    
    return res
end
# Sometimes we have additional data corruption that leads to a non-matching trajectory lengths
# In order to identify the spurious data, we compare the central energies for every trajectory
# We then remove duplicates  and/or non-matching central energies from the replica with the longer trajectory
function find_first_duplicated_central_energies(S0,traj_lengths)
    Smin,Smax = extrema(filter(isfinite, S0))
    nreplicas = first(size(S0))
    S0_theory = collect(range(Smin,Smax,length=nreplicas))
    traj_min  = minimum(traj_lengths)
    for i in 1:traj_min
        if !isapprox(S0_theory,sort(S0[:,i]))
            # identify duplicate indices and then
            # remove those that correspond to the longer trajectories
            # ans make sure that only duplicates have occured in pairs
            inds = approx_duplicate_indices(S0[:,i])
            return i, inds
        end
    end
end
function remove_central_energy_trajectories!(S0,traj_step,replica_indices)
    for replica_ind in replica_indices
        S0[replica_ind,traj_step:end-1] = S0[replica_ind,traj_step+1:end]
        # The last entry is then set to infinity, since it corresponds to missing data
        # and we choose 'Inf' to represent that because it autoimatically places the entry
        # at the end of the array when sorting it.
        S0[replica_ind,end] = Inf 
    end
    return S0
end
function remove_non_matching_trajectories_in_replicas(fid_repeat)
    S0 = read_non_matching_trajectory(fid_repeat,key="S0")
    remove_non_matching_trajectories_in_replicas(S0)
end
function remove_non_matching_trajectories_in_replicas(S0::Matrix{T}) where T
    traj_lengths = dropdims(count(isfinite, S0, dims=2),dims=2)

    ans = find_first_duplicated_central_energies(S0, traj_lengths)
    while !isnothing(ans)
        traj_step, indices = ans
        rm_indices = filter(i -> traj_lengths[i] > minimum(traj_lengths),  indices)
        S0 = remove_central_energy_trajectories!(S0,traj_step,rm_indices)
        traj_lengths = dropdims(count(isfinite, S0, dims=2),dims=2)
        ans = find_first_duplicated_central_energies(S0, traj_lengths)
    end

    # In the end we check that all trajectories up to minimum(traj_lengths) are matching the 
    # initial central energies and the remaining entries all correspond to missing data which 
    # we encode as 'Inf'
    min_traj, max_traj = extrema(traj_lengths)
    Smin,Smax = extrema(filter(isfinite, S0))
    nreplicas = first(size(S0))
    S0_theory = collect(range(Smin,Smax,length=nreplicas))
    for traj in 1:min_traj
        @assert isapprox(S0_theory,sort(S0[:,traj]))
    end
    for traj in min_traj+1:max_traj
        @assert all(isequal(Inf),S0[:,traj])
    end
    # Then we only return the good, matching trajectories
    return S0[:,1:min_traj]
end