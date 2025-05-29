using Pkg; Pkg.activate(".")
using LLRParsing
using HiRepOutputCleaner
using DelimitedFiles

# specified directory
path          = "/home/fabian/Documents/Physics/Data/" 
metadata      = "metadata/runs.csv"
# temporary directories after cleaning
path          = "tmp" 
metadata      = "metadata/tmp.csv"
h5file        = "data_assets/test_Nt5.hdf5"
h5file_sorted = "data_assets/test_Nt5_sorted.hdf5"

function llr_alldirs_hdf5(path,metadata_file,h5file)
    metadata = readdlm(metadata_file,',',skipstart=1)
    for (run,clean) in eachrow(metadata)
        dir = joinpath(path,run)
        if clean
            tmpdir = "./tmp/$(basename(dir))/"
            ispath(tmpdir) || mkpath(tmpdir)
            clean_llr_directory(dir,tmpdir;checkpoint_pattern=nothing,last_ranges=nothing,warn=false)
            #llr_dir_hdf5(tmpdir,h5file)
            #rm(tmpdir,recursive=true)
        else
            llr_dir_hdf5(dir,h5file)
        end
    end
end

isfile(h5file) && rm(h5file)
isfile(h5file_sorted) && rm(h5file_sorted)
llr_alldirs_hdf5(path,metadata,h5file)
sort_by_central_energy_to_hdf5(h5file, h5file_sorted)