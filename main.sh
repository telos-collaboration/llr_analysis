path="."
sdir="scripts/done" 
julia scripts/instantiate.jl

# h5 files
Nt4_hdf5="data_assets/Sp4_Nt4_sorted.hdf5"
Nt5_hdf5="data_assets/Sp4_Nt5_sorted.hdf5"
Nt6_hdf5="data_assets/Sp4_Nt6_sorted.hdf5"

# parsing
julia --project="." $sdir/parse_llr.jl --path $path --metadata "metadata/runsNt4.csv" --h5file_unsorted "data_assets/Sp4_Nt4.hdf5" --h5file $Nt4_hdf5 --clean false
julia --project="." $sdir/parse_llr.jl --path $path --metadata "metadata/runsNt5.csv" --h5file_unsorted "data_assets/Sp4_Nt5.hdf5" --h5file $Nt5_hdf5 --clean false
julia --project="." $sdir/parse_llr.jl --path $path --metadata "metadata/runsNt6.csv" --h5file_unsorted "data_assets/Sp4_Nt6.hdf5" --h5file $Nt6_hdf5 --clean false

# analysis
julia --project="." $sdir/trajectory_overview.jl --h5file $Nt4_hdf5 --plot_dir "data_assets/plots/overview"
julia --project="." $sdir/trajectory_overview.jl --h5file $Nt5_hdf5 --plot_dir "data_assets/plots/overview"
julia --project="." $sdir/trajectory_overview.jl --h5file $Nt6_hdf5 --plot_dir "data_assets/plots/overview"

julia --project="." $sdir/an_histories.jl --h5file $Nt4_hdf5 --plot_dir "data_assets/plots/an_trajectories"
julia --project="." $sdir/an_histories.jl --h5file $Nt5_hdf5 --plot_dir "data_assets/plots/an_trajectories"
julia --project="." $sdir/an_histories.jl --h5file $Nt6_hdf5 --plot_dir "data_assets/plots/an_trajectories"

julia --project="." $sdir/free_energy.jl --h5file $Nt4_hdf5 --plot_dir "assets/plots/free_energy"
julia --project="." $sdir/free_energy.jl --h5file $Nt5_hdf5 --plot_dir "assets/plots/free_energy"

julia --project="." $sdir/compare_volumes.jl --h5file $Nt4_hdf5 --plot_file "assets/plots/an_volume_Nt4.pdf" --title "\$a_n\$ volume dependence for \$N_t=4\$" --xmin 0.5695 --xmax 0.573  --ymin 7.337 --ymax 7.343
julia --project="." $sdir/compare_volumes.jl --h5file $Nt5_hdf5 --plot_file "assets/plots/an_volume_Nt5.pdf" --title "\$a_n\$ volume dependence for \$N_t=5\$" --xmin 0.5889 --xmax 0.5898 --ymin 7.488 --ymax 7.491
julia --project="." $sdir/compare_volumes.jl --h5file $Nt6_hdf5 --plot_file "assets/plots/an_volume_Nt6.pdf" --title "\$a_n\$ volume dependence for \$N_t=6\$" --xmin 0.6015 --xmax 0.6025 --ymin 7.615 --ymax 7.623

julia --project="." $sdir/double_gaussian_fit.jl --h5file $Nt4_hdf5 --plot_dir "assets/plots/plaquette_distribution" --xmin 0.567  --xmax 0.577  --beta_min 7.337 --beta_max 7.343
julia --project="." $sdir/double_gaussian_fit.jl --h5file $Nt5_hdf5 --plot_dir "assets/plots/plaquette_distribution" --xmin 0.5885 --xmax 0.5905 --beta_min 7.488 --beta_max 7.492

julia --project="." $sdir/double_gaussian_volumes.jl --h5file $Nt4_hdf5 --plotfile "assets/plots/plaquette_distribution_Nt4_volumes.pdf" --xmin 0.567  --xmax 0.577  --beta_min 7.337 --beta_max 7.343
julia --project="." $sdir/double_gaussian_volumes.jl --h5file $Nt5_hdf5 --plotfile "assets/plots/plaquette_distribution_Nt5_volumes.pdf" --xmin 0.5885 --xmax 0.5905 --beta_min 7.488 --beta_max 7.492

julia --project="." $sdir/critical_beta.jl --h5file $Nt4_hdf5 --outfile "data_assets/critical_beta_Nt4.csv"     --beta_min 7.337 --beta_max 7.343
julia --project="." $sdir/critical_beta.jl --h5file $Nt4_hdf5 --outfile "data_assets/critical_beta_Nt4_2:1.csv" --beta_min 7.337 --beta_max 7.343 --relative_peak_height 2
julia --project="." $sdir/critical_beta.jl --h5file $Nt5_hdf5 --outfile "data_assets/critical_beta_Nt5.csv"     --beta_min 7.488 --beta_max 7.492
julia --project="." $sdir/critical_beta.jl --h5file $Nt5_hdf5 --outfile "data_assets/critical_beta_Nt5_2:1.csv" --beta_min 7.488 --beta_max 7.492 --relative_peak_height 2

julia --project="." $sdir/plot_beta.jl --plotfile "assets/plots/critical_beta_volumes_Nt4.pdf" "data_assets/critical_beta_Nt4.csv" "data_assets/critical_beta_Nt4_2:1.csv"
julia --project="." $sdir/plot_beta.jl --plotfile "assets/plots/critical_beta_volumes_Nt5.pdf" "data_assets/critical_beta_Nt5.csv" "data_assets/critical_beta_Nt5_2:1.csv"

julia --project="." $sdir/tables.jl --h5file $Nt5_hdf5 --outfile "assets/tables/runs.tex"

