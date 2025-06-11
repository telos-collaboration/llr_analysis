using Pkg; Pkg.activate(".")
using LLRParsing
using HiRepOutputCleaner
using DelimitedFiles

function llr_alldirs_hdf5(path,metadata_file,h5file,cleaned=true)
    metadata = readdlm(metadata_file,',',skipstart=1)
    tmpdir   = "./tmp/"
    if !cleaned
        for run in eachrow(String.(metadata))
            dir = joinpath(path,run[1])
            p   = joinpath(tmpdir,basename(dir))
            ispath(p) || mkpath(p)
            clean_llr_directory(dir,p;checkpoint_pattern=nothing,last_ranges=nothing,warn=false)
        end
    end
    for dir in readdir(tmpdir,join=true)
        llr_dir_hdf5(dir,h5file)
    end
end

# Nt=5
path          = "/home/fabian/Documents/Physics/Data/" 
metadata      = "metadata/runs.csv"
h5file        = "data_assets/Sp4_Nt5.hdf5"
h5file_sorted = "data_assets/Sp4_Nt5_sorted.hdf5"

isfile(h5file) && rm(h5file)
isfile(h5file_sorted) && rm(h5file_sorted)
llr_alldirs_hdf5(path,metadata,h5file)
sort_by_central_energy_to_hdf5(h5file, h5file_sorted)