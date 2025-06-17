path="."
julia scripts/instantiate.jl
# Nt6
julia scripts/done/parse_llr.jl --path $path --metadata "metadata/runsNt6.csv" --h5file_unsorted "data_assets/Sp4_Nt6.hdf5" --h5file "data_assets/Sp4_Nt6_sorted.hdf5" --clean false
julia scripts/done/compare_volumes.jl --h5file "data_assets/Sp4_Nt6_sorted.hdf5" --plot_file "assets/an_volume_Nt6.pdf" --title "\$a_n\$ volume dependence for \$N_t=6\$" --xmin 0.6015 --xmax 0.6025 --ymin 7.615 --ymax 7.623