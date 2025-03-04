# This function uses the Parsers.jl package to quickly parse 
# a string that contains 'n' doubles seperated by spaces.
# Then result is stored in pre-allocated array 'tmp'
function _parse_float_string!(array,string;n)
    opts = Parsers.Options(delim=' ', ignorerepeated=true)
    io = IOBuffer(string)
    for i in 1:n
        array[i] = Parsers.parse(Float64, io, opts)
    end
end
# Parse temporal Polyakov loop from the logfile
# (i.e. we only parse direction '0')
# I extract the full complex loop
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
function importance_sampling_file_to_hdf5(file_in,h5file)
    Nt, Nl, plaq, beta, poly = parse_importance_sampling(file_in)
    ensemble = "ImportanceSampling/$(Nt)x$(Nl)/$beta"

    h5write(h5file,joinpath(ensemble,"plaquette"),plaq)
    h5write(h5file,joinpath(ensemble,"polyakov_loop"),poly)
    h5write(h5file,joinpath(ensemble,"Nt"),Nt)
    h5write(h5file,joinpath(ensemble,"Nl"),Nl)
    h5write(h5file,joinpath(ensemble,"beta"),beta)
end
function importance_sampling_dir_hdf5(basedir,h5file;outname="output_file")
    isfile(h5file) && rm(h5file)
    # Iterate recursively through all directories in basedir
    # Check if there is a matching output file
    # If yes: Save data to the hdf5 file 
    for (root,dirs,files) in walkdir(basedir)
        if outname ∈ files
            infile = joinpath(root,outname)
            importance_sampling_file_to_hdf5(infile,h5file)
        end
    end
end

