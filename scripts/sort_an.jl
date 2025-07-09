using LLRParsing
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
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
    h5file = args["h5file_unsorted"]
    h5file_sorted = args["h5file"]
    ispath(dirname(h5file_sorted)) || mkpath(dirname(h5file_sorted))
    isfile(h5file_sorted) && rm(h5file_sorted)

    return sort_by_central_energy_to_hdf5(h5file, h5file_sorted)
end
main()
