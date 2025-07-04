sdir="scripts"
julia scripts/instantiate.jl

# h5 files
Nt4_hdf5="data_assets/Sp4_Nt4_sorted.hdf5"
Nt5_hdf5="data_assets/Sp4_Nt5_sorted.hdf5"
Nt6_hdf5="data_assets/Sp4_Nt6_sorted.hdf5"

# parsing
#julia --project="." $sdir/parse_llr.jl --metadata "metadata/runs.csv" --h5file_unsorted "tmp/Sp4_Nt4.hdf5" --Nt 4
#julia --project="." $sdir/parse_llr.jl --metadata "metadata/runs.csv" --h5file_unsorted "tmp/Sp4_Nt5.hdf5" --Nt 5
#julia --project="." $sdir/parse_llr.jl --metadata "metadata/runs.csv" --h5file_unsorted "tmp/Sp4_Nt6.hdf5" --Nt 6

# sorting of a_n
#julia --project="." $sdir/sort_an.jl --h5file_unsorted "tmp/Sp4_Nt4.hdf5" --h5file $Nt4_hdf5
#julia --project="." $sdir/sort_an.jl --h5file_unsorted "tmp/Sp4_Nt5.hdf5" --h5file $Nt5_hdf5
#julia --project="." $sdir/sort_an.jl --h5file_unsorted "tmp/Sp4_Nt6.hdf5" --h5file $Nt6_hdf5

## analysis
#julia --project="." $sdir/trajectory_overview.jl --h5file $Nt4_hdf5 --plot_dir "data_assets/plots/overview"
#julia --project="." $sdir/trajectory_overview.jl --h5file $Nt5_hdf5 --plot_dir "data_assets/plots/overview"
#julia --project="." $sdir/trajectory_overview.jl --h5file $Nt6_hdf5 --plot_dir "data_assets/plots/overview"
##julia --project="." $sdir/one_trajectory_overview.jl --h5file $Nt5_hdf5 --run_name "5x56_128replicas" --plot_file "data_assets/plots/overview/5x56_128replicas.pdf"

#julia --project="." $sdir/an_histories.jl --h5file $Nt4_hdf5 --plot_dir "data_assets/plots/an_trajectories"
#julia --project="." $sdir/an_histories.jl --h5file $Nt5_hdf5 --plot_dir "data_assets/plots/an_trajectories"
#julia --project="." $sdir/an_histories.jl --h5file $Nt6_hdf5 --plot_dir "data_assets/plots/an_trajectories"
#julia --project="." $sdir/one_an_history.jl --h5file $Nt6_hdf5 --plotfile "data_assets/plots/an_trajectories/$run" --run_name $run

#julia --project="." $sdir/free_energy.jl --h5file $Nt4_hdf5 --plot_dir "assets/plots/free_energy"
#julia --project="." $sdir/free_energy.jl --h5file $Nt5_hdf5 --plot_dir "assets/plots/free_energy"
#
#julia --project="." $sdir/entropy.jl --h5file $Nt4_hdf5 --plot_file "assets/plots/entropy_Nt4.pdf"
#julia --project="." $sdir/entropy.jl --h5file $Nt5_hdf5 --plot_file "assets/plots/entropy_Nt5.pdf"
#
#julia --project="." $sdir/compare_volumes.jl --h5file $Nt4_hdf5 --plot_file "assets/plots/an_volume_Nt4.pdf" --title "\$a_n\$ volume dependence for \$N_t=4\$"
#julia --project="." $sdir/compare_volumes.jl --h5file $Nt5_hdf5 --plot_file "assets/plots/an_volume_Nt5.pdf" --title "\$a_n\$ volume dependence for \$N_t=5\$"
#julia --project="." $sdir/compare_volumes.jl --h5file $Nt6_hdf5 --plot_file "assets/plots/an_volume_Nt6.pdf" --title "\$a_n\$ volume dependence for \$N_t=6\$"
#
##julia --project="." $sdir/tables.jl --h5file $Nt5_hdf5 --outfile "assets/tables/runs.tex"
#
#julia --project="." $sdir/critical_beta.jl --h5file $Nt4_hdf5 --outfile "data_assets/critical_beta_Nt4.csv"     
#julia --project="." $sdir/critical_beta.jl --h5file $Nt4_hdf5 --outfile "data_assets/critical_beta_Nt4_2:1.csv"  --relative_peak_height 2
#julia --project="." $sdir/critical_beta.jl --h5file $Nt5_hdf5 --outfile "data_assets/critical_beta_Nt5.csv"     
#julia --project="." $sdir/critical_beta.jl --h5file $Nt5_hdf5 --outfile "data_assets/critical_beta_Nt5_2:1.csv"  --relative_peak_height 2
#
#julia --project="." $sdir/double_gaussian_fit.jl --h5file $Nt4_hdf5 --plot_dir "assets/plots/plaquette_distribution"
#julia --project="." $sdir/double_gaussian_fit.jl --h5file $Nt5_hdf5 --plot_dir "assets/plots/plaquette_distribution"
#
#julia --project="." $sdir/double_gaussian_volumes.jl --h5file $Nt4_hdf5 --plotfile "assets/plots/plaquette_distribution_Nt4_volumes.pdf"
#julia --project="." $sdir/double_gaussian_volumes.jl --h5file $Nt5_hdf5 --plotfile "assets/plots/plaquette_distribution_Nt5_volumes.pdf"
#
#julia --project="." $sdir/plot_beta.jl --plotfile "assets/plots/critical_beta_volumes_Nt4.pdf" "data_assets/critical_beta_Nt4.csv" "data_assets/critical_beta_Nt4_2:1.csv"
#julia --project="." $sdir/plot_beta.jl --plotfile "assets/plots/critical_beta_volumes_Nt5.pdf" "data_assets/critical_beta_Nt5.csv" "data_assets/critical_beta_Nt5_2:1.csv"