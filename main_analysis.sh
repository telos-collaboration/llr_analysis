h5file="data_assets/Sp4_Nt5_sorted.hdf5"
plotdir="data_assets/overview_Nt5"

julia scripts/trajectory_overview.jl --h5file $h5file --plot_dir $plotdir
