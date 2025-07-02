using LLRParsing
using DelimitedFiles
using ArgParse

function parse_skip(str)
    spl = split(chop(str, head=1),',')
    return all(isempty,spl) ? String[] : String.(spl)
end
function parse_full(dir,skip, h5file)
    s = parse_skip.(skip)
    llr_dir_hdf5(dir,h5file;skip_repeats=s)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--metadata"
            help = "CSV file containing the metadata of runs to be analysed"
            required = true
        "--h5file_unsorted"
            help = "Where to write the resulting HDF5 file containing the unsorted results"
            required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    metadata_file = args["metadata"] 
    h5file        = args["h5file_unsorted"]
    metadata      = readdlm(metadata_file,',',String,skipstart=1)

    ispath(dirname(h5file)) || mkpath(dirname(h5file))
    isfile(h5file) && rm(h5file)

    for row in eachrow(metadata)
        run, s = row
        parse_full(run,s, h5file)
    end
end
main()