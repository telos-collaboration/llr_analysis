using Pkg; Pkg.activate(".")
using LLRParsing
using HiRepOutputCleaner
using DelimitedFiles
using ArgParse

function clean_all_llr_dirs(path,metadata_file;target="./tmp/")
    metadata = readdlm(metadata_file,',',skipstart=1)
    for run in eachrow(String.(metadata))
        dir = joinpath(path,run[1])
        p   = joinpath(target,basename(dir))
        ispath(p) || mkpath(p)
        clean_llr_directory(dir,p;checkpoint_pattern=nothing,last_ranges=nothing,warn=false)
    end
end
function llr_alldirs_hdf5(path,metadata_file,h5file;clean=false,tmpdir="./tmp/")
    metadata = readdlm(metadata_file,',',skipstart=1)
    runs     = joinpath.(path,String.(metadata))
    if clean
        clean_all_llr_dirs(path,metadata_file;target=tmpdir)
        runs = readdir(tmpdir,join=true)
    end
    for dir in runs
        llr_dir_hdf5(dir,h5file)
    end
    if clean
        rm(tmpdir,recursive=true)
    end
end
function parse_full(path,metadata, h5file, h5file_sorted,clean)
    isfile(h5file) && rm(h5file)
    isfile(h5file_sorted) && rm(h5file_sorted)
    llr_alldirs_hdf5(path,metadata,h5file; clean)
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
        "--clean"
            help = "Run cleaning script to remove unfinished runs from logs"
            arg_type = Bool
            default = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    path          = args["path"] 
    metadata      = args["metadata"] 
    h5file        = args["h5file_unsorted"]
    h5file_sorted = args["h5file"] 
    clean         = args["clean"] 
    parse_full(path,metadata, h5file, h5file_sorted,clean)
end
main()