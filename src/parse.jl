function _parse_float_string!(array,string;n) 
    opts = Parsers.Options(delim=' ', ignorerepeated=true)
    io = IOBuffer(string)
    for i in 1:n
        array[i] = Parsers.parse(Float64, io, opts)
    end
end
function polyakov_loop(file)
    polyakov_loop = ComplexF64[]
    tmp = zeros(ComplexF32,2)
    pattern = "[FUND_POLYAKOV][0]Polyakov direction 0 = "
    i = length(pattern)
    for line in eachline(file)
        if startswith(line,pattern)
            _parse_float_string!(tmp,line[i:end];n=2)
            append!(polyakov_loop,tmp[1] + im*tmp[2])
        end
    end
    return polyakov_loop
end
function parse_beta(file)
    for line in eachline(file)
        if startswith(line,"[INIT][0]beta=")
            beta = parse(Float64,line[length("[INIT][0]beta=")+1:end])
            return beta
        end
    end
end
function parse_importance_sampling(file)
    Nt, Nl = latticesize(file)[1:2]
    plaq = plaquettes(file)
    beta = parse_beta(file)
    poly = polyakov_loop(file)
    return Nt, Nl, plaq, beta, poly 
end
function importance_sampling_hdf5(file_in,h5file)
    Nt, Nl, plaq, beta, poly = parse_importance_sampling(file_in)
    ensemble = "ImportanceSampling/$(Nt)x$(Nl)/$beta"

    h5write(h5file,joinpath(ensemble,"plaquette"),plaq)
    h5write(h5file,joinpath(ensemble,"polyakov_loop"),poly)
    h5write(h5file,joinpath(ensemble,"Nt"),Nt)
    h5write(h5file,joinpath(ensemble,"Nl"),Nl)
    h5write(h5file,joinpath(ensemble,"beta"),beta)
end