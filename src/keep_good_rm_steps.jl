using HDF5
using LLRParsing
# This function determines all trajectory numbers that appear in the output file 
# Note, that the number itself can be misleading: The starting number is typically parsed from the ouput file of the preceeding run.
# If the preceeding run has errored then the configuration will still be named after the last number appearing in the out_0 file even if the configuration has never been written to disk
function traj_numbers(file)
    pattern = r"^\[MAIN\]\[0\](Newton Raphson|Robbins Monro) sequence \#([0-9]+): generated"
    numbers = Int[]
    types   = String[]
    for line in eachline(file)
        if occursin(pattern,line)
            m = match(pattern,line)
            push!(types,m.captures[1])
            push!(numbers,parse(Int,m.captures[2]))
        end
    end
    return numbers,types
end
# The trajectory number in the output file cannot be trusted. To circumvent this issue I determine all trajectory numbers that indicate that a configuration has been 
# read from or written to the disk. With this I can identify the Robbins-Monro trajectory numbers that should be included in the analysis
function traj_read_saved(file)
    NR_or_RM = r"(NR |)Trajectory"
    pattern  = r"Configuration.*n([0-9]+)\].*(read|saved)"
    type     = ""
    numbers    = Int[]
    types      = String[]
    read_saved = String[]
    for line in eachline(file)
        if occursin(NR_or_RM,line)
            m = match(NR_or_RM,line)
            t = m.captures[1]
            type = t == "" ? "RM" : "NR" 
        end
        if occursin(pattern,line)
            m = match(pattern,line)
            number        = m.captures[1]
            read_or_saved = m.captures[2]
            push!(numbers,parse(Int,number)) 
            push!(types,type) 
            push!(read_saved,read_or_saved)
        end
    end
    # Once we are done, we need to take care of the mislabelling, when we change from NR to RM
    # When reading in the new configuration the loop above did not xte know that a RM step is coming up
    for i in eachindex(read_saved, types, numbers)
        i == 1 && continue
        read_NR  = read_saved[i-1] == "read"  && types[i-1] == "NR"
        saved_RM = read_saved[i]   == "saved" && types[i]   == "RM"
        if read_NR && saved_RM 
            types[i-1] = "RM"
        end
    end
    return read_saved, types, numbers
end
# Identify healthy trajectory numbers in any given run 
# Note, that David's script currently rename the configuration before reading by decreasing the configuration number by 1
# The first occurence of a (read|save) pattern should be a save pattern
# We declare a Robbins Monro trajectry number good if it appears between a 'read' and a 'save' pattern 
function find_good_RM_trajectory_numbers(read_saved, types, numbers)
    good_trajectories = Int[]
    # if we haven't done any NR steps, the first time we save to disk will be a RM iteration without
    # There will not be a matching configuration that we have read from disk
    if read_saved[1] == "saved" && types[1] == "RM" 
        r = 0:numbers[1]
        append!(good_trajectories,r)
    end
    # Otherwise look for RM pairs of first reading and then saving
    for i in eachindex(read_saved, types, numbers)
        i == 1 && continue
        if types[i] == "RM"
            if read_saved[i] == "saved" && read_saved[i-1] == "read" 
                r = numbers[i-1]+1:numbers[i]-1
                append!(good_trajectories,r)
            end
        end
    end
    return good_trajectories
end
function find_good_RM_trajectory_numbers(file)
    read_saved, types, rs_number = traj_read_saved(file)
    good_trajectories = find_good_RM_trajectory_numbers(read_saved, types, rs_number)
    return good_trajectories
end
function good_rm_indices(trajectories,types,good_trajectories)
    good_indices = Int[]
    for i in eachindex(trajectories)
        if types[i] == "Newton Raphson"
            push!(good_indices,i)
        elseif types[i] == "Robbins Monro" && trajectories[i] ∈ good_trajectories
            push!(good_indices,i)
        end
    end
    return good_indices
end
function parse_llr_corrupted(file)
    dS0, S0, plaq, a, is_rm, S0_fxa, a_fxa, poly = parse_llr(file)
    trajectories, types = traj_numbers(file)
    if !(size(S0) == size(plaq) == size(a) == size(is_rm))
        @show file
        @show size(S0), size(plaq), size(a), size(is_rm)
        @assert size(S0) == size(plaq) == size(a) == size(is_rm)
    end
    if !(size(S0) == size(trajectories) == size(types))
        @show size(S0), size(trajectories), size(types)
        @assert size(S0) == size(trajectories) == size(types)
    end
    good_trajectories = find_good_RM_trajectory_numbers(file) 
    inds = good_rm_indices(trajectories,types,good_trajectories)
    return dS0, S0[inds], plaq[inds], a[inds], is_rm[inds], S0_fxa, a_fxa, poly
end
