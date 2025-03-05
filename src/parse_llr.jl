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
function parse_S0(file)
    pattern = "[SWAP][10]New Rep Par S0 = "
    pos1    = length(pattern)
    S0      = Float64[]
    is_fxa = false
    for line in eachline(file)
        if startswith(line,"[SYSTEM][0]Process finalized.")
            is_fxa = false
        end
        if startswith(line,"[MAIN][0]Robins Monro update done.")
            is_fxa = true
        end        
        if !is_fxa && startswith(line,pattern)
            pos2 = first(findnext("dS",line,pos1))
            append!(S0,parse(Float64,line[pos1:pos2-1]))
        end
    end
    return S0
end
function parse_llr_plaquette(file)
    pattern = r"^\[MAIN\]\[0\](NR )*Plaq a fixed ([0-9]+.[0-9]+)"
    plaq    = Float64[]
    is_rm   = Bool[]
    for line in eachline(file)
        if occursin(pattern,line)
            m   = match(pattern,line)
            str = m.captures
            append!(plaq,parse(Float64,str[2]))
            append!(is_rm,isnothing(str[1]))
        end
    end
    return plaq, is_rm
end
function parse_a(file)
    pattern = r"^\[MAIN\]\[0\](NR )*<a_rho\(.+\)>= ([0-9]+.[0-9]+)"
    a       = Float64[]
    is_rm   = Bool[]
    for line in eachline(file)
        if occursin(pattern,line)
            m   = match(pattern,line)
            str = m.captures
            append!(a,parse(Float64,str[2]))
            append!(is_rm,isnothing(str[1]))
        end
    end
    return a, is_rm
end
function parse_fixeda_S0_a_dS(file;S0_last,a_last,dS_last)
    is_fxa = false
    S0 = Float64[S0_last]
    a  = Float64[a_last]
    rx = r" S0 ([0-9]+.[0-9]+),  a  ([0-9]+.[0-9]+) , dS ([0-9]+.[0-9]+)"
    for line in eachline(file)
        if startswith(line,"[SYSTEM][0]Process finalized.")
            is_fxa = false
        end
        if startswith(line,"[MAIN][0]Robins Monro update done.")
            is_fxa = true
        end
        if is_fxa && startswith(line,"[llr:setreplica][0]New LLR Param:")
            vals = match(rx,line).captures 
            append!(S0,parse(Float64,vals[1]))
            append!(a ,parse(Float64,vals[2]))
            @assert parse(Float64,vals[3]) == dS_last
        end
    end
    # remove the last entry, because no measurement has been performed after the last fixed-a update
    return S0[1:end-1], a[1:end-1]
end
function _parse_data!(array,string;n)
    opts = Parsers.Options(delim=' ', ignorerepeated=true)
    io = IOBuffer(string)
    for i in 1:n
        array[i] = Parsers.parse(Float64, io, opts)
    end
end
function parse_fun_polyakov_loop(file)
    poly    = ComplexF64[]
    pattern = "[FUND_POLYAKOV][0]Polyakov direction 0 = "
    pos     = length(pattern)
    vals    = zeros(2)
    for line in eachline(file)
        if startswith(line,pattern)
            _parse_data!(vals,line[pos:end];n=2)
            append!(poly, vals[1] + im*vals[2])
        end
    end
    return poly
end
function parse_llr_no_fxa(file)
    dS0 = parse_dS0(file)
    S0  = parse_S0(file)
    plaq, is_rmP = parse_llr_plaquette(file)
    a, is_rma   = parse_a(file) 
    @assert is_rmP == is_rma
    return dS0, S0, plaq, a, is_rma
end
function parse_llr_full(file)
    dS0, S0, plaq, a, is_rm = parse_llr_no_fxa(file)
    S0_fxa, a_fxa = parse_fixeda_S0_a_dS(file;S0_last=S0[end],a_last=a[end],dS_last=dS0)
    poly = parse_fun_polyakov_loop(file)
    return dS0, S0, plaq, a, is_rm, S0_fxa, a_fxa, poly
end
function parse_llr_quick(file)
    """ This is equvivalent to parse_llr_full(file) but performs only one pass over the file  """
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