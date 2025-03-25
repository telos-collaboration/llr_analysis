using HDF5
using LLRParsing
# This function determines all trajectory numbers that appear in the output file 
# Note, that the number itself can be misleading: The starting number is typically parsed from the ouput file of the preceeding run.
# If the preceeding run has errored then the configuration will still be named after the last number appearing in the out_0 file even if the configuration has never been written to disk
function traj_numbers(file)
    pattern = r"^\[MAIN\]\[0\](Newton Raphson|Robbins Monro) sequence \#([0-9]+): generated"
    numbers = String[]
    for line in eachline(file)
        if occursin(pattern,line)
            m = match(pattern,line)
            step = m.captures[1]*" "*m.captures[2]
            push!(numbers,step)
        end
    end
    return numbers
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
function find_good_RM_trajectory_numbers(read_saved)
    good_trajectories = Int[]
    old_type = ""
    for id in read_saved
        m = match(r"([0-9]+) (read|saved) (RM|NR)",id)
        number = parse(Int,m.captures[1]) 
        rs     = m.captures[2]
        RM_NR  = m.captures[3]
        @show number, rs, RM_NR
    end
end

# Showcase to files where a remnant data corruption occurs
# for i in [0,2]
#     file = "output/LLRout/LLR_5x72_152/0/Rep_$i/out_0"
#     dS0, S0, plaq, a, is_rm, S0_fxa, a_fxa, poly = parse_llr(file)
#     numbers = traj_numbers(file)
#     read_saved = traj_read_saved(file)
#     @show i, size(S0), size(a), size(plaq), size(is_rm)
# end

file = "output/LLRout/LLR_5x72_152/0/Rep_0/out_0"
file = "output/LLRout/LLR_5x72_152/0/Rep_2/out_0"
dS0, S0, plaq, a, is_rm, S0_fxa, a_fxa, poly = parse_llr(file)
numbers = traj_numbers(file)
read_saved, types, rs_number = traj_read_saved(file)
@show size(read_saved), size(types), size(rs_number)
