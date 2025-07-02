using LLRParsing
using DelimitedFiles
using ArgParse

function parse_skip(str)
    spl = split(chop(str, head=1),',')
    return all(isempty,spl) ? String[] : String.(spl)
end
function llr_alldirs_hdf5(path,metadata_file,h5file;tmpdir="./tmp/")
    ispath(dirname(h5file)) || mkpath(dirname(h5file))
    metadata = readdlm(metadata_file,',',String,skipstart=1)
    dirs = metadata[:,1]
    skip = parse_skip.(metadata[:,2])
    runs = joinpath.(path,dirs)
    for (dir,s) in zip(runs,skip)
        llr_dir_hdf5(dir,h5file;skip_repeats=s)
    end
end
function parse_full(path,metadata, h5file, h5file_sorted)
    isfile(h5file) && rm(h5file)
    isfile(h5file_sorted) && rm(h5file_sorted)
    llr_alldirs_hdf5(path,metadata,h5file)
    sort_by_central_energy_to_hdf5(h5file, h5file_sorted)
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--metadata"
            help = "CSV file containing the metadata of runs to be analysed"
            required = true
        "--path"
            help = "Directory relative to file locations specified in metadata file"
            required = true
        "--h5file"
            help = "Where to write the resulting HDF5 file containing the sorted results"
            required = true
        "--h5file_unsorted"
            help = "Where to write the resulting HDF5 file containing the unsorted results"
            required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    path          = args["path"] 
    metadata      = args["metadata"] 
    h5file        = args["h5file_unsorted"]
    h5file_sorted = args["h5file"] 
    parse_full(path,metadata, h5file, h5file_sorted)
end
main()