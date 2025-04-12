using LLRParsing
using HDF5

function read_initial_a(h5dset,run;repeat)
    replicas = keys(h5dset[run][repeat])
    Nl = h5dset[run]["Nl"][]
    Nt = h5dset[run]["Nt"][]
    a0 = [ read(h5dset[run][repeat][repl],"a0") for repl in replicas]
    S0 = [ h5dset[run][repeat][repl]["S0_sorted"][1] for repl in replicas]
    u0 = S0/(6*Nl^3*Nt)
    @show Nl, Nt
    return a0, u0
end

dir   = "/home/fabian/Documents/Physics/Data/DataSunbird/Sunbird/LLR/LLR_4x20_40" 
dir   = "/home/fabian/Dokumente/Physics/Data/DataSunbird/Sunbird/LLR/LLR_4x20_40" 
file  = "test_small.hdf5"
file0 = "test_small_sorted.hdf5"
#llr_dir_hdf5(dir,file)
#sort_by_central_energy_to_hdf5(file, file0)

h5dset = h5open(file0)
runs   = keys(h5dset)
a0, u0 = read_initial_a(h5dset,runs[1];repeat="0")

using Plots
plt = scatter(u0,a0,label="intial a")

a, S   = LLRParsing.a_vs_central_action_repeats(h5dset,runs[1];ind=nothing)[1:2]
LLRParsing.a_vs_central_action_plot!(plt,a[:,1],zero(a[:,1]),S[:,1],4,20;lens=false)
LLRParsing.a_vs_central_action_plot!(plt,a[:,2],zero(a[:,2]),S[:,2],4,20;lens=false)

a, S   = LLRParsing.a_vs_central_action_repeats(h5dset,runs[1];ind=1)[1:2]
LLRParsing.a_vs_central_action_plot!(plt,a[:,1],zero(a[:,1]),S[:,1],4,20;lens=false)
LLRParsing.a_vs_central_action_plot!(plt,a[:,2],zero(a[:,2]),S[:,2],4,20;lens=false)