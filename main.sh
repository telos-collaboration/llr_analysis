path="."
julia scripts/instantiate.jl

# parsing
julia scripts/done/parse_llr.jl --path $path --metadata "metadata/runsNt5.csv" --h5file_unsorted "data_assets/Sp4_Nt5.hdf5" --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --clean false
julia scripts/done/parse_llr.jl --path $path --metadata "metadata/runsNt4.csv" --h5file_unsorted "data_assets/Sp4_Nt4.hdf5" --h5file "data_assets/Sp4_Nt4_sorted.hdf5" --clean false

# analysis
julia scripts/done/trajectory_overview.jl --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --plot_dir "assets/overview_Nt5"
julia scripts/done/trajectory_overview.jl --h5file "data_assets/Sp4_Nt4_sorted.hdf5" --plot_dir "assets/overview_Nt4"
julia scripts/done/compare_volumes.jl --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --plot_file "assets/an_volume_Nt5.pdf" --title "\$a_n\$ volume dependence for \$N_t=5\$" --xmin 0.5889 --xmax 0.5898 --ymin 7.488 --ymax 7.491
julia scripts/done/compare_volumes.jl --h5file "data_assets/Sp4_Nt4_sorted.hdf5" --plot_file "assets/an_volume_Nt4.pdf" --title "\$a_n\$ volume dependence for \$N_t=4\$" --xmin 0.5695 --xmax 0.573  --ymin 7.337 --ymax 7.343
julia scripts/done/double_gaussian_fit.jl --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --plot_dir "assets/plaquette_distribution_Nt5" --xmin 0.5885 --xmax 0.5905 --beta_min 7.488 --beta_max 7.492
julia scripts/done/double_gaussian_fit.jl --h5file "data_assets/Sp4_Nt4_sorted.hdf5" --plot_dir "assets/plaquette_distribution_Nt4" --xmin 0.567  --xmax 0.577  --beta_min 7.337 --beta_max 7.343
julia scripts/done/tables.jl --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --outfile "assets/runs.tex"
julia scripts/done/double_gaussian_volumes.jl --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --plotfile "assets/plaquette_distribution_Nt5_volumes.pdf" --xmin 0.5885 --xmax 0.5905 --beta_min 7.488 --beta_max 7.492
julia scripts/done/double_gaussian_volumes.jl --h5file "data_assets/Sp4_Nt4_sorted.hdf5" --plotfile "assets/plaquette_distribution_Nt4_volumes.pdf" --xmin 0.567  --xmax 0.577  --beta_min 7.337 --beta_max 7.343
julia scripts/done/critical_beta.jl --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --outfile "data_assets/critical_beta_Nt5.csv" --beta_min 7.488 --beta_max 7.492
julia scripts/done/critical_beta.jl --h5file "data_assets/Sp4_Nt4_sorted.hdf5" --outfile "data_assets/critical_beta_Nt4.csv" --beta_min 7.337 --beta_max 7.343