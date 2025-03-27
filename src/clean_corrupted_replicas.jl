# Utilities for 
# 1) Finding all indices for unique elements in an array v  
# 2) Finding the values of all items of v that occur more than once
# 3) Finding the inidces of all items of v that occur more than once allowing for slight deviations
unique_indices(v) = unique(i -> v[i], eachindex(v))
duplicate_vals(v) = v[setdiff(eachindex(v),unique_indices(v))]
approx_duplicate_indices(v) = findall(i-> any(isapprox(v[i]),duplicate_vals(v)), eachindex(v))

function read_non_matching_trajectory(fid,::Type{T};key) where T
    replicas     = keys(fid)
    traj_lengths = [length(fid["$rep/$key"]) for rep in replicas]
    nmin, nmax   = extrema(traj_lengths)
    res          = fill(typemax(T),(length(replicas),nmax))
    for (i,rep) in enumerate(replicas)
        val = read(fid,"$rep/$key")
        res[i,1:length(val)] = val
    end
    return res
end
# Sometimes we have additional data corruption that leads to a non-matching trajectory lengths
# In order to identify the spurious data, we compare the central energies for every trajectory
# We then remove duplicates  and/or non-matching central energies from the replica with the longer trajectory
function find_first_duplicated_central_energies(S0,traj_lengths;i_min=1)
    Smin,Smax = extrema(filter(isfinite, S0))
    nreplicas = first(size(S0))
    S0_theory = collect(range(Smin,Smax,length=nreplicas))
    traj_min  = minimum(traj_lengths)
    for i in i_min:traj_min
        if !isapprox(S0_theory,sort(S0[:,i]))
            # identify duplicate indices and then
            # remove those that correspond to the longer trajectories
            # ans make sure that only duplicates have occured in pairs
            inds = approx_duplicate_indices(S0[:,i])
            return i, inds
        end
    end
end
function remove_trajectories!(data::Matrix{T},traj_step,replica_indices) where T
    for replica_ind in replica_indices
        data[replica_ind,traj_step:end-1] = data[replica_ind,traj_step+1:end]
        # The last entry is then set to infinity (the typemax for Float64), since it corresponds to missing data
        # and we choose typemax(T) ('Inf' for Float64) to represent that, because it automatically places the entry
        # at the end of the array when sorting it.
        data[replica_ind,end] = typemax(T)
    end
end
function remove_all_trajectories!(data::Matrix{T},traj_step_array,replica_indices_array) where T
    for (traj_step,replica_indices) in zip(traj_step_array,replica_indices_array)
        remove_trajectories!(data,traj_step,replica_indices)
    end
end
function remove_non_matching_trajectories_in_replicas(fid_repeat)
    S0 = read_non_matching_trajectory(fid_repeat,Float64,key="S0")
    remove_non_matching_trajectories_in_replicas(S0)
end
function remove_non_matching_trajectories_in_replicas(S0::Matrix{T}) where T
    traj_lengths = dropdims(count(isfinite, S0, dims=2),dims=2)

    # keep track of trajectory numbers and indices of non-matching trajectories
    # (such that they can later on be applied to the other observables)
    steps = Int[]
    inds  = Array{Int}[]

    ans = find_first_duplicated_central_energies(S0, traj_lengths)
    while !isnothing(ans)
        traj_step, indices = ans
        push!(steps,traj_step)
        push!(inds ,indices)
        rm_indices = filter(i -> traj_lengths[i] > minimum(traj_lengths),  indices)
        remove_trajectories!(S0,traj_step,rm_indices)
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
    return S0[:,1:min_traj], steps, inds
end