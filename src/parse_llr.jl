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
        dir_dict[repeat] = replica_dirs 
    end
    return repeat_dirs, dir_dict
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
function parse_llr(file)
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

    for line in eachline(file)
    
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
    repeats, replica_dirs = get_repeat_and_replica_dirs(dir)
    files      = _all_files_from_dict(dir,replica_dirs)
    N_repeats  = length(repeats)
    # assure the global lattice parameters are identical for all repeats and replicas
    N_replicas = only(unique([length(replica_dirs[r]) for r in repeats]))
    Nt = only(unique(first.(latticesize.(files))))  
    Nl = only(unique(last.(latticesize.(files))))

    name = "$(Nt)x$(Nl)_$(N_repeats)repeats_$(N_replicas)replicas"*suffix
    write(fid,joinpath(name,"N_repeats"),N_repeats)
    write(fid,joinpath(name,"N_replicas"),N_replicas)
    write(fid,joinpath(name,"Nt"),Nt)
    write(fid,joinpath(name,"Nl"),Nl)

    @showprogress desc="parsing $name" for repeat in repeats
        for rep in replica_dirs[repeat]
            file = joinpath(dir, repeat,rep,"out_0")
            dS0, S0, plaq, a, is_rm, S0_fxa, a_fxa, poly = parse_llr(file)
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
function sort_by_central_energy_to_hdf5(h5file_in,h5file_out)
    h5dset = h5open(h5file_in)
    for run in keys(h5dset)
        run == "6x72_26repeats_48replicas" && continue
        N_replicas = read(h5dset[run],"N_replicas")
        N_repeats  = read(h5dset[run],"N_repeats")
        # read all last elements for a and the central action
        for j in 1:N_repeats
            ntraj = length(h5dset[run]["$(j-1)/Rep_0/is_rm"])
            a = zeros(N_replicas,ntraj)
            S = zeros(N_replicas,ntraj)
            p = zeros(N_replicas,ntraj)
            for i in 1:N_replicas
                a[i,:] = h5dset[run]["$(j-1)/Rep_$(i-1)/a"][] 
                S[i,:] = h5dset[run]["$(j-1)/Rep_$(i-1)/S0"][]
                p[i,:] = h5dset[run]["$(j-1)/Rep_$(i-1)/plaq"][]
            end
            ## Sort by the central action to account for different swaps
            for j in 1:ntraj
                perm = sortperm(S[:,j])
                S[:,j] = S[perm,j]
                a[:,j] = a[perm,j]
                p[:,j] = p[perm,j]
            end
            ## make sure that the sorted central action alwas matches
            for i in 1:N_replicas
                @assert allequal(S[i,:])
                h5write(h5file_out,joinpath(run,"$(j-1)","Rep_$(i-1)","S0_sorted"),  S[i,:])
                h5write(h5file_out,joinpath(run,"$(j-1)","Rep_$(i-1)","a_sorted"),   a[i,:])
                h5write(h5file_out,joinpath(run,"$(j-1)","Rep_$(i-1)","plaq_sorted"),p[i,:])
                h5write(h5file_out,joinpath(run,"$(j-1)","Rep_$(i-1)","dS0"),h5read(h5file_in,joinpath(run,"$(j-1)","Rep_$(i-1)","dS0")))
                h5write(h5file_out,joinpath(run,"$(j-1)","Rep_$(i-1)","is_rm"),h5read(h5file_in,joinpath(run,"$(j-1)","Rep_$(i-1)","is_rm")))
            end
        end
        h5write(h5file_out,joinpath(run,"N_replicas"),N_replicas)
        h5write(h5file_out,joinpath(run,"N_repeats"),N_repeats)
        h5write(h5file_out,joinpath(run,"Nt"),h5read(h5file_in,joinpath(run,"Nt")))
        h5write(h5file_out,joinpath(run,"Nl"),h5read(h5file_in,joinpath(run,"Nl")))
    end
end