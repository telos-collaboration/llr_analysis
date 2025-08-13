# Utilities for
# 1) Finding all indices for unique elements in an array v
# 2) Finding the values of all items of v that occur more than once
# 3) Finding the inidces of all items of v that occur more than once allowing for slight deviations
unique_indices(v) = unique(i -> v[i], eachindex(v))
duplicate_vals(v) = v[setdiff(eachindex(v), unique_indices(v))]
approx_duplicate_indices(v) =
    findall(i -> any(isapprox(v[i]), duplicate_vals(v)), eachindex(v))

function read_non_matching_trajectory(fid, ::Type{T}; key) where {T}
    replicas = keys(fid)
    traj_lengths = [length(fid["$rep/$key"]) for rep in replicas]
    nmin, nmax = extrema(traj_lengths)
    res = fill(typemax(T), (length(replicas), nmax))
    for (i, rep) in enumerate(replicas)
        val = read(fid, "$rep/$key")
        res[i, 1:length(val)] = val
    end
    return res
end
# Sometimes we have additional data corruption that leads to a non-matching trajectory lengths
# In order to identify the spurious data, we compare the central energies for every trajectory
# We then remove duplicates  and/or non-matching central energies from the replica with the longer trajectory
function find_first_duplicated_central_energies(S0, traj_lengths; i_min = 1)
    Smin, Smax = extrema(filter(isfinite, S0))
    nreplicas = first(size(S0))
    S0_theory = collect(range(Smin, Smax, length = nreplicas))
    traj_min = minimum(traj_lengths)
    for i in i_min:traj_min
        if !isapprox(S0_theory, sort(S0[:, i]))
            # identify duplicate indices and then
            # remove those that correspond to the longer trajectories
            # ans make sure that only duplicates have occured in pairs
            inds = approx_duplicate_indices(S0[:, i])
            return i, inds
        end
    end
    return nothing, nothing
end
