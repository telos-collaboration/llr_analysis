using LLRParsing
using DelimitedFiles
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--outfile"
        help = "Where to save the combined file"
        required = true
        "files"
        help = "CSV file(s) to be combines"
        required = true
        nargs = '+'
    end
    return parse_args(s)
end
function combine_csv(outfile, files)
    csvs = readdlm.(files, ',', comments = true, header = true)
    data = getindex.(csvs, 1)
    headers = getindex.(csvs, 2)
    @assert allequal(headers)
    new_csv = vcat(first(headers), data...)

    io = open(outfile, "w")
    print_provenance_csv(io)
    writedlm(io, new_csv, ',')
    return close(io)
end
function main()
    args = parse_commandline()
    return combine_csv(args["outfile"], args["files"])
end
main()
