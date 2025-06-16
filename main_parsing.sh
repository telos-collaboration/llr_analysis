# Nt=5
julia scripts/parse_llr.jl --path "/home/fabian/Documents/Physics/Data/" --metadata "metadata/runs.csv" --h5file_unsorted "data_assets/Sp4_Nt5.hdf5" --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --clean true
# Nt=4
julia scripts/parse_llr.jl --path "/home/fabian/Documents/Physics/Data/" --metadata "metadata/runsNt4.csv" --h5file_unsorted "data_assets/Sp4_Nt4.hdf5" --h5file "data_assets/Sp4_Nt4_sorted.hdf5" --clean false