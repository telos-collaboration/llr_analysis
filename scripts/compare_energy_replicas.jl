using HDF5
# Utilities for 
# 1) Finding all indices for unique elements in an array v  
# 2) Finding the values of all items of v that occur more than once
# 3) Finding the inidces of all items of v that occur more than once
nique_indices(v) = unique(i -> v[i], eachindex(v))
duplicate_vals(v) = v[setdiff(eachindex(v),unique_indices(v))]
duplicate_indices(v) = findall(i-> v[i] ∈ duplicate_vals(v), eachindex(v))

function all_central_enrgies(fid)
    replicas     = keys(fid)
    traj_lengths = [length(fid["$rep/S0"]) for rep in replicas]
    nmin, nmax   = extrema(traj_lengths)
    S0 = zeros((length(replicas),nmax)) .+ Inf
    for (i,rep) in enumerate(replicas)
        val = read(fid_repeat,"$rep/S0")
        S0[i,1:length(val)] = val
    end
    return S0, replicas, traj_lengths
end
# Sometimes we have additional data corruption that leads to a non-matching trajectory lengths
# In order to identify the spurious data, we compare the central energies for every trajectory
# We then remove duplicates  and/or non-matching central energies from the replica with the longer trajectory
function find_first_duplicated_central_energies(S0,replicas,traj_lengths)
    Smin,Smax = extrema(filter(isfinite, S0))
    S0_theory = collect(range(Smin,Smax,length=length(replicas)))
    traj_min  = minimum(traj_lengths)
    for i in 1:traj_min
        if !isapprox(S0_theory,sort(S0[:,i]))
            # identify duplicate indices and then
            # remove those that correspond to the longer trajectories
            # ans make sure that only duplicates have occured in pairs
            inds = duplicate_indices(S0[:,i])
            return i, inds
        end
    end
end
function remove_central_energy_trajectories!(S0,traj_step,replica_indices)
    for replica_ind in replica_indices
        S0[replica_ind,traj_step:end-1] = S0[replica_ind,traj_step+1:end]
    end
    return S0
end

fid = h5open("output/test.hdf5")
fid_repeat = fid["5x64_19repeats_95replicas/13/"]

S0, replicas, traj_lengths = all_central_enrgies(fid_repeat)
traj_step, indices = find_first_duplicated_central_energies(S0, replicas, traj_lengths)
S0 = remove_central_energy_trajectories!(S0,traj_step,indices)
traj_step, indices = find_first_duplicated_central_energies(S0, replicas, traj_lengths)

Smin,Smax = extrema(filter(isfinite, S0))
S0_theory = collect(range(Smin,Smax,length=length(replicas)))
