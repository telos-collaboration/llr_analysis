using Pkg; Pkg.activate(".")
using LLRParsing
using HiRepOutputCleaner

path = "/home/fabian/Documents/Physics/Data/" 
dirs = [
    "DataSunbird/Sunbird/LLR/LLR_5x48_48_v2",
    "DataTursa/LLR/LLR_5x80_64",
    "DataCSD/CSD3/LLR_5x64_95",
    "DataCSD/CSD3/LLR_5x72_95",
]

h5file        = "output/test_Nt5.hdf5"
h5file_sorted = "output/test_Nt5_sorted.hdf5"

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