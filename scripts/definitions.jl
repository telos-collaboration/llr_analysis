using HDF5
using LLRParsing
using ArgParse
function read_all_llr_therm(fid, runs)
    llr_therm = Int[]
    llr_meas = Int[]
    for r in runs
        N_replicas = read(fid[r], "N_replicas")
        repeats = read(fid[r], "repeats")
        for repeat in repeats, rep in 0:(N_replicas - 1)
            append!(llr_therm, read(fid, joinpath(r, repeat, "Rep_$rep", "llr_therm")))
            append!(llr_meas, read(fid, joinpath(r, repeat, "Rep_$rep", "llr_meas")))
        end
    end
    # assert that we always were using the same parameters
    llr_therm = only(unique(llr_therm))
    llr_meas = only(unique(llr_meas))
    return llr_therm, llr_meas
end
function definitions(h5file_Nt4, h5file_Nt5, outfile)
    fid_Nt4 = h5open(h5file_Nt4)
    fid_Nt5 = h5open(h5file_Nt5)

    runs_Nt4 = filter(!startswith("provenance"), keys(fid_Nt4))
    runs_Nt5 = filter(!startswith("provenance"), keys(fid_Nt5))

    ratio_Nt4 = [ read(fid_Nt4, joinpath(r, "Nl")) / read(fid_Nt4, joinpath(r, "Nt")) for r in runs_Nt4 ]
    ratio_Nt5 = [ read(fid_Nt5, joinpath(r, "Nl")) / read(fid_Nt5, joinpath(r, "Nt")) for r in runs_Nt5 ]

    max_aspect_ratio_Nt4 = Int(maximum(ratio_Nt4))
    max_aspect_ratio_Nt5 = Int(maximum(ratio_Nt5))
    llr_therm, llr_meas = read_all_llr_therm(fid_Nt5, runs_Nt5)

    ispath(dirname(outfile)) || mkpath(dirname(outfile))
    @show outfile
    io = open(outfile, "w")
    print_provenance_tex(io)
    write(io, "\\newcommand\\MaxAspectRatioNtFour{$max_aspect_ratio_Nt4}\n")
    write(io, "\\newcommand\\MaxAspectRatioNtFive{$max_aspect_ratio_Nt5}\n")
    write(io, "\\newcommand\\llrtherm{$llr_therm}\n")
    write(io, "\\newcommand\\llrmeas{$llr_meas}\n")
    return close(io)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file_Nt4"
        help = "HDF5 file containing the unsorted results at Nt=4"
        required = true
        "--h5file_Nt5"
        help = "HDF5 file containing the unsorted results at Nt=5"
        required = true
        "--outfile"
        help = "TeX file for writing the definitions to disk"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    h5file_Nt4 = args["h5file_Nt4"]
    h5file_Nt5 = args["h5file_Nt5"]
    outfile = args["outfile"]
    return definitions(h5file_Nt4, h5file_Nt5, outfile)
end
main()
