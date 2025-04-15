using LLRParsing
using HDF5
using Plots
gr(fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

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

dir1  = "/home/fabian/Documents/Physics/Data/DataSunbird/Sunbird/LLR/LLR_4x20_40" 
dir2  = "/home/fabian/Documents/Physics/Data/DataSunbird/Sunbird/LLR/LLR_4x20_40_v1" 
file  = "test_small.hdf5"
file0 = "test_small_sorted.hdf5"
isfile(file)  && rm(file)
isfile(file0) && rm(file0)
llr_dir_hdf5(dir1,file)
llr_dir_hdf5(dir2,file)
sort_by_central_energy_to_hdf5(file, file0)
h5dset = h5open(file0)
runs   = keys(h5dset)


a0, u0 = read_initial_a(h5dset,runs[2];repeat="0")
plt1   = scatter(u0,a0,label="intial a")
a, S   = LLRParsing.a_vs_central_action_repeats(h5dset,runs[2];ind=1)[1:2]
LLRParsing.a_vs_central_action_plot!(plt1,a[:,1],zero(a[:,1]),S[:,1],4,20;lens=false)
a, S   = LLRParsing.a_vs_central_action_repeats(h5dset,runs[2];ind=10)[1:2]
LLRParsing.a_vs_central_action_plot!(plt1,a[:,1],zero(a[:,1]),S[:,1],4,20;lens=false)

a0, u0 = read_initial_a(h5dset,runs[2];repeat="1")
plt2   = scatter(u0,a0,label="intial a")
a, S   = LLRParsing.a_vs_central_action_repeats(h5dset,runs[2];ind=1)[1:2]
LLRParsing.a_vs_central_action_plot!(plt2,a[:,2],zero(a[:,2]),S[:,2],4,20;lens=false)
a, S   = LLRParsing.a_vs_central_action_repeats(h5dset,runs[2];ind=10)[1:2]
LLRParsing.a_vs_central_action_plot!(plt2,a[:,2],zero(a[:,2]),S[:,2],4,20;lens=false)

a0, u0 = read_initial_a(h5dset,runs[1];repeat="0")
plt3   = scatter(u0,a0,label="intial a")
a, S   = LLRParsing.a_vs_central_action_repeats(h5dset,runs[1];ind=1)[1:2]
LLRParsing.a_vs_central_action_plot!(plt3,a[:,1],zero(a[:,1]),S[:,1],4,20;lens=false)
a, S   = LLRParsing.a_vs_central_action_repeats(h5dset,runs[1];ind=10)[1:2]
LLRParsing.a_vs_central_action_plot!(plt3,a[:,1],zero(a[:,1]),S[:,1],4,20;lens=false)

xl = xlims(plt1)
yl = ylims(plt1)
for plt in [plt1,plt2,plt3]
    plot!(plt,xlims=xl,ylims=yl)
end

plot!(plt1,title="n_or=0, domains=1, 1 NR step (orange) vs 10 NR steps (green)")
plot!(plt2,title="n_or=4, domains=1, 1 NR step (orange) vs 10 NR steps (green)")
plot!(plt3,title="n_or=0, domains=4, 1 NR step (orange) vs 10 NR steps (green)")
savefig(plt1,"no_or_no_dd.pdf")
savefig(plt2,"only_or.pdf")
savefig(plt3,"only_dd.pdf")
display(plt1)
display(plt2)
display(plt3)