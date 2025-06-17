path="."
julia scripts/instantiate.jl

julia scripts/done/parse_llr.jl --path $path --metadata "metadata/runsNt5.csv" --h5file_unsorted "data_assets/Sp4_Nt5.hdf5" --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --clean false
julia scripts/done/parse_llr.jl --path $path --metadata "metadata/runsNt4.csv" --h5file_unsorted "data_assets/Sp4_Nt4.hdf5" --h5file "data_assets/Sp4_Nt4_sorted.hdf5" --clean false