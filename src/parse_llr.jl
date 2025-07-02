function get_repeat_and_replica_dirs(base_dir,skip_repeats=String[])
    # Obtain all directories containing repeats and sort them
    repeat_dirs = filter( str -> all(isdigit, str), readdir(base_dir))
    sort!(repeat_dirs,lt=natural)
    if !isempty(skip_repeats)
        repeat_dirs = filter(i-> i ∉ skip_repeats,repeat_dirs)
    end    
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
function parse_initial_a(file)
    a0 = NaN
    pattern = "[MAIN][0]LLR Initial a"
    for line in eachline(file)
        if startswith(line,pattern)
            a0 = parse(Float64,line[length(pattern)+1:end])
            return a0
        end
    end
    return a0
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
                if !isempty(S0) && !isempty(S0) 
                    is_fxa = true
                    append!(S0_fxa,S0[end])
                    append!(a_fxa,a[end])
                end
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
function llr_dir_hdf5(dir,h5file;suffix="",skip_repeats=String[]) 
    fid = h5open(h5file,"cw")

    # get all repeats and replicas and store that information for future use
    replica_dirs = get_repeat_and_replica_dirs(dir,skip_repeats)
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
            a0   = parse_initial_a(file)
            dS0, S0, plaq, a, is_rm, S0_fxa, a_fxa, poly = parse_llr(file)
            write(fid,joinpath(name,repeat,rep,"dS0"),dS0)
            write(fid,joinpath(name,repeat,rep,"S0"),S0)
            write(fid,joinpath(name,repeat,rep,"a0"),a0)
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
    h5dset = h5open(h5file_in,"r")
    runs   = keys(h5dset)
    close(h5dset)
    for run in runs
        sort_by_central_energy_to_hdf5_run(h5file_in,h5file_out,run)
    end
end
function sort_by_central_energy_to_hdf5_run(h5file_in,h5file_out,run)
    h5dset     = h5open(h5file_in,"r")
    h5dset_out = h5open(h5file_out,"cw")

    N_replicas = read(h5dset[run],"N_replicas")
    N_repeats  = read(h5dset[run],"N_repeats")
    repeats    = read(h5dset[run],"repeats")        
    # read all last elements for a and the central action
    for j in repeats
        # read data for all replicas
        a     = read_non_matching_trajectory(h5dset[run][j],Float64;key="a")
        p     = read_non_matching_trajectory(h5dset[run][j],Float64;key="plaq")
        is_rm = read_non_matching_trajectory(h5dset[run][j],Bool   ;key="is_rm")
        S     = read_non_matching_trajectory(h5dset[run][j],Float64;key="S0")
        # Check if we have any mismatches of the unsorted central energies
        data_healthy = all(allequal ,eachslice(sort(S,dims=1),dims=1))
        if !data_healthy
            traj_lengths = dropdims(count(isfinite, S, dims=2),dims=2)
            last_healthy_traj_p1, inds = find_first_duplicated_central_energies(S, traj_lengths)
            S = S[:,1:last_healthy_traj_p1-1]
            @warn "Run $run, repeat $j: Discarded data after step $(last_healthy_traj_p1-1)"
            data_healthy = all(allequal ,eachslice(sort(S,dims=1),dims=1))
        end

        ntraj = dropdims(count(isfinite, S, dims=2),dims=2)
        n_traj_min, n_traj_max = extrema(ntraj)
        ## Sort by the central action to account for different swaps
        for j in 1:n_traj_min
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
        # make sure that the sorted central action alwas matches, if not, discard the repeat
        @assert data_healthy
        for i in 1:N_replicas
            dset    = create_group(h5dset_out, joinpath(run,"$j","Rep_$(i-1)"))
            dset_in = h5dset[joinpath(run,"$j","Rep_$(i-1)")]
            dS0 = read(dset_in,"dS0")
            a0  = read(dset_in,"a0")
            write(dset,"S0_sorted",  S[i,:])
            write(dset,"a_sorted",   a[i,:])
            write(dset,"plaq_sorted",p[i,:])
            write(dset,"is_rm"      ,is_rm[i,:])
            write(dset,"dS0"        ,dS0)
            write(dset,"a0"         ,a0)
        end
    end
    write(h5dset_out,joinpath(run,"N_replicas"),N_replicas)
    write(h5dset_out,joinpath(run,"N_repeats"),N_repeats)
    write(h5dset_out,joinpath(run,"repeats"),repeats)
    write(h5dset_out,joinpath(run,"Nt"),h5read(h5file_in,joinpath(run,"Nt")))
    write(h5dset_out,joinpath(run,"Nl"),h5read(h5file_in,joinpath(run,"Nl")))

    close(h5dset)
    close(h5dset_out)
end