function get_repeat_and_replica_dirs(base_dir)
    # Obtain all directories containing repeats and sort them
    repeat_dirs = filter( str -> all(isdigit, str), readdir(base_dir))
    sort!(repeat_dirs,lt=natural)
    
    dir_dict = Dict{String,Vector{String}}()
    for repeat in repeat_dirs
        rx = r"Rep_[0-9]+"
        repeat_path  = joinpath(base_dir,repeat)  
        replica_dirs = filter(startswith(rx),readdir(repeat_path))
        sort!(replica_dirs,lt=natural)
        # check if there exists an ouput file for every replica in this repeat
        # If not, then skip this repeat. This scenario rarely happens and only
        # has been seen in thermalisation
        files = joinpath.(Ref(base_dir),Ref(repeat),replica_dirs,Ref("out_0"))
        if all(isfile, files)
            dir_dict[repeat] = replica_dirs
        else
            @warn "directory $(basename(base_dir)), repeat $repeat: some output files are missing/empty"
        end
    end
    return dir_dict
end
function _all_files_from_dict(dir,replica_dirs)
    files = AbstractString[]
    for repeat in keys(replica_dirs), rep in replica_dirs[repeat]
        push!(files, joinpath(dir,repeat,rep,"out_0"))
    end
    return files
end
function parse_dS0(file)
    dS0 = NaN
    pattern = "[MAIN][0]LLR Delta S"
    for line in eachline(file)
        if startswith(line,pattern)
            dS0 = parse(Float64,line[length(pattern)+1:end])
            return dS0
        end
    end
    return dS0
end
function _parse_data!(array,string;n)
    opts = Parsers.Options(delim=' ', ignorerepeated=true)
    io = IOBuffer(string)
    for i in 1:n
        array[i] = Parsers.parse(Float64, io, opts)
    end
end
function parse_llr(file; skiplines=Int[])
    pattern_poly = "[FUND_POLYAKOV][0]Polyakov direction 0 = "
    patternS0 = "[SWAP][10]New Rep Par S0 = "
    patternPl = r"^\[MAIN\]\[0\](NR )*Plaq a fixed ([0-9]+.[0-9]+)"
    pattern_a = r"^\[MAIN\]\[0\](NR )*<a_rho\(.+\)>= ([0-9]+.[0-9]+)"
    pos_poly  = length(pattern_poly)
    posS0 = length(patternS0)
    dS0 = parse_dS0(file)
    rx = r" S0 ([0-9]+.[0-9]+),  a  ([0-9]+.[0-9]+) , dS ([0-9]+.[0-9]+)"
    
    is_rm  = Bool[]
    plaq   = Float64[]
    S0     = Float64[]
    a      = Float64[]
    poly   = ComplexF64[]
    
    S0_fxa = Float64[]
    a_fxa  = Float64[]
    
    tmp_poly = zeros(2)
    is_fxa = false

    # keep track of line number so that we can skip them if specified by skiplines
    for (line_no,line) in enumerate(eachline(file))            
        # check if we want to skip the current line
        line_no ∈ skiplines && continue
        # then continue with the usual parsing
        if startswith(line,"[SYSTEM][0]Process finalized.")
            is_fxa = false
        end
        if startswith(line,"[MAIN][0]")
            if startswith(line,"[MAIN][0]Robins Monro update done.")
                is_fxa = true
                append!(S0_fxa,S0[end])
                append!(a_fxa,a[end])
            end
            if occursin(patternPl,line)
                m   = match(patternPl,line)
                str = m.captures
                append!(plaq,parse(Float64,str[2]))
                append!(is_rm,isnothing(str[1]))
            end
            if occursin(pattern_a,line)
                m   = match(pattern_a,line)
                str = m.captures
                append!(a,parse(Float64,str[2]))
            end
        end
        if !is_fxa && startswith(line,patternS0)
            pos2 = first(findnext("dS",line,posS0))
            append!(S0,parse(Float64,line[posS0:pos2-1]))
        end
        if is_fxa && startswith(line,"[llr:setreplica][0]New LLR Param:")
            vals = match(rx,line).captures 
            append!(S0_fxa,parse(Float64,vals[1]))
            append!(a_fxa ,parse(Float64,vals[2]))
            @assert parse(Float64,vals[3]) == dS0
        end
        if startswith(line,pattern_poly)
            _parse_data!(tmp_poly,line[pos_poly:end];n=2)
            append!(poly, tmp_poly[1] + im*tmp_poly[2])
        end
    end
    return dS0, S0, plaq, a, is_rm, S0_fxa[1:end-1], a_fxa[1:end-1], poly
end
function llr_dir_hdf5(dir,h5file;suffix="") 
    fid = h5open(h5file,"cw")

    # get all repeats and replicas and store that information for future use
    replica_dirs = get_repeat_and_replica_dirs(dir)
    repeats      = sort(collect(keys(replica_dirs)),lt=natural)
    files        = _all_files_from_dict(dir,replica_dirs)
    N_repeats    = length(repeats)
    if isempty(repeats)
        @warn "No non-emtpy logfiles available for $dir"
        return
    end
    # assure the global lattice parameters are identical for all repeats and replicas
    N_replicas = only(unique([length(replica_dirs[r]) for r in repeats]))
    Nt = only(unique(first.(latticesize.(files))))  
    Nl = only(unique(last.(latticesize.(files))))

    name = "$(Nt)x$(Nl)_$(N_repeats)repeats_$(N_replicas)replicas"*suffix
    write(fid,joinpath(name,"N_repeats"),N_repeats)
    write(fid,joinpath(name,"N_replicas"),N_replicas)
    write(fid,joinpath(name,"repeats"),repeats)
    write(fid,joinpath(name,"Nt"),Nt)
    write(fid,joinpath(name,"Nl"),Nl)

    @showprogress desc="parsing $name" for repeat in repeats
        for rep in replica_dirs[repeat]
            file = joinpath(dir, repeat,rep,"out_0")
            dS0, S0, plaq, a, is_rm, S0_fxa, a_fxa, poly = parse_llr_corrupted(file)
            write(fid,joinpath(name,repeat,rep,"dS0"),dS0)
            write(fid,joinpath(name,repeat,rep,"S0"),S0)
            write(fid,joinpath(name,repeat,rep,"plaq"),plaq)
            write(fid,joinpath(name,repeat,rep,"a"),a)
            write(fid,joinpath(name,repeat,rep,"is_rm"),is_rm)
            write(fid,joinpath(name,repeat,rep,"S0_fxa"),S0_fxa)
            write(fid,joinpath(name,repeat,rep,"a_fxa"),a_fxa)
            write(fid,joinpath(name,repeat,rep,"poly"),poly)
        end
    end
    close(fid)
end
function sort_by_central_energy_to_hdf5(h5file_in,h5file_out;skip_ens=nothing)
    h5dset     = h5open(h5file_in,"r")
    h5dset_out = h5open(h5file_out,"cw")
    @showprogress desc="sorting $h5file_in" for run in keys(h5dset)
        if !isnothing(skip_ens) 
            run ∈ skip_ens && continue
        end
        N_replicas = read(h5dset[run],"N_replicas")
        N_repeats  = read(h5dset[run],"N_repeats")
        repeats    = read(h5dset[run],"repeats")        
        # read all last elements for a and the central action
        for j in repeats
            ntraj1 = [ length(h5dset[run]["$j/Rep_$i/a"])     for i in 0:N_replicas-1]
            ntraj2 = [ length(h5dset[run]["$j/Rep_$i/S0"])    for i in 0:N_replicas-1]
            ntraj3 = [ length(h5dset[run]["$j/Rep_$i/plaq"])  for i in 0:N_replicas-1]
            ntraj4 = [ length(h5dset[run]["$j/Rep_$i/is_rm"]) for i in 0:N_replicas-1]
            @assert ntraj1 == ntraj2 == ntraj3 == ntraj4            
            n_traj_min, n_traj_max = extrema(ntraj1)
            # I need to deal with those later
            # I want to remove extra trajectories that correspond to non matching replicas
            if n_traj_min < n_traj_max
                @warn "Run $run, repeat $j: Non-matching trajectory length across replicas. Non-matching entries will be discarded."
            end
            # read data for all replicas
            a     = read_non_matching_trajectory(h5dset[run][j],Float64;key="a")
            p     = read_non_matching_trajectory(h5dset[run][j],Float64;key="plaq")
            is_rm = read_non_matching_trajectory(h5dset[run][j],Bool   ;key="is_rm")
            S     = read_non_matching_trajectory(h5dset[run][j],Float64;key="S0")
            ## Sort by the central action to account for different swaps
            for j in 1:n_traj_max
                perm = sortperm(S[:,j])
                S[:,j] = S[perm,j]
                a[:,j] = a[perm,j]
                p[:,j] = p[perm,j]
                is_rm[:,j] = is_rm[perm,j]
            end
            a = a[:,1:n_traj_min]
            p = p[:,1:n_traj_min]
            S = S[:,1:n_traj_min]
            is_rm = is_rm[:,1:n_traj_min]
            ## make sure that the sorted central action alwas matches, if not, discard the repeat
            keep_repeat = true
            for i in 1:N_replicas 
                keep_repeat = keep_repeat*allequal(S[i,:])
            end
            if keep_repeat
                for i in 1:N_replicas
                    dset    = create_group(h5dset_out, joinpath(run,"$j","Rep_$(i-1)"))
                    dset_in = h5dset[joinpath(run,"$j","Rep_$(i-1)")]
                    dS0 = read(dset_in,"dS0")
                    write(dset,"S0_sorted",  S[i,:])
                    write(dset,"a_sorted",   a[i,:])
                    write(dset,"plaq_sorted",p[i,:])
                    write(dset,"is_rm"      ,is_rm[i,:])
                    write(dset,"dS0"        ,dS0)
                end
            else
                @warn "Run $run, repeat $j: Central energies do not always match => discarded "
            end
        end
        write(h5dset_out,joinpath(run,"N_replicas"),N_replicas)
        write(h5dset_out,joinpath(run,"N_repeats"),N_repeats)
        write(h5dset_out,joinpath(run,"repeats"),repeats)
        write(h5dset_out,joinpath(run,"Nt"),h5read(h5file_in,joinpath(run,"Nt")))
        write(h5dset_out,joinpath(run,"Nl"),h5read(h5file_in,joinpath(run,"Nl")))
    end
end