using Pkg; Pkg.activate(".")
using LLRParsing
using HiRepOutputCleaner
using DelimitedFiles

path = "/home/fabian/Documents/Physics/Data/" 
dirs = readdlm("metadata/runs.csv",',',skipstart=1)

h5file        = "data_assets/test_Nt5.hdf5"
h5file_sorted = "data_assets/test_Nt5_sorted.hdf5"

function llr_alldirs_hdf5(path,dirs,h5file)
    for dir in joinpath.(Ref(path),dirs)
        tmpdir = "./tmp/$(basename(dir))/"
        ispath(tmpdir) || mkpath(tmpdir)
        clean_llr_directory(dir,tmpdir;checkpoint_pattern=nothing,last_ranges=nothing,warn=false)
        llr_dir_hdf5(tmpdir,h5file)
        rm(tmpdir,recursive=true)
    end
end

isfile(h5file) && rm(h5file)
isfile(h5file_sorted) && rm(h5file_sorted)
llr_alldirs_hdf5(path,dirs,h5file)
sort_by_central_energy_to_hdf5(h5file, h5file_sorted)