julia scripts/trajectory_overview.jl --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --plot_dir "data_assets/overview_Nt5"
julia scripts/trajectory_overview.jl --h5file "data_assets/Sp4_Nt4_sorted.hdf5" --plot_dir "data_assets/overview_Nt4"

julia scripts/compare_volumes.jl --h5file "data_assets/Sp4_Nt5_sorted.hdf5" --plot_file "assets/an_volume_Nt5.pdf" --title "\$a_n\$ volume dependence for \$N_t=5\$" --xmin 0.5889 --xmax 0.5898 --ymin 7.488 --ymax 7.491
julia scripts/compare_volumes.jl --h5file "data_assets/Sp4_Nt4_sorted.hdf5" --plot_file "assets/an_volume_Nt4.pdf" --title "\$a_n\$ volume dependence for \$N_t=4\$" --xmin 0.5695 --xmax 0.573  --ymin 7.337 --ymax 7.343