using Pkg; Pkg.activate(".")
using LLRParsing
using HiRepOutputCleaner
using DelimitedFiles

function clean_all_llr_dirs(path,metadata_file;target="./tmp/")
    metadata = readdlm(metadata_file,',',skipstart=1)
    for run in eachrow(String.(metadata))
        dir = joinpath(path,run[1])
        p   = joinpath(target,basename(dir))
        ispath(p) || mkpath(p)
        clean_llr_directory(dir,p;checkpoint_pattern=nothing,last_ranges=nothing,warn=false)
    end
end
function llr_alldirs_hdf5(path,metadata_file,h5file;clean=false)
    metadata = readdlm(metadata_file,',',skipstart=1)
    runs     = joinpath.(path,String.(metadata))
    if clean
        tmpdir = "./tmp/"
        clean_all_llr_dirs(path,metadata_file;target=tmpdir)
        runs = readdir(tmpdir,join=true)
    end
    for dir in runs
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
llr_alldirs_hdf5(path,metadata,h5file; clean=true)
sort_by_central_energy_to_hdf5(h5file, h5file_sorted)

# Nt=4
#path          = "/home/fabian/Documents/Physics/Data/" 
#metadata      = "metadata/runsNt4.csv"
#h5file        = "data_assets/Sp4_Nt4.hdf5"
#h5file_sorted = "data_assets/Sp4_Nt4_sorted.hdf5"
#
#isfile(h5file) && rm(h5file)
#isfile(h5file_sorted) && rm(h5file_sorted)
#llr_alldirs_hdf5(path,metadata,h5file; clean=false)
#sort_by_central_energy_to_hdf5(h5file, h5file_sorted)