using HDF5
using Plots
using Statistics

file = "output/SU3_llr_david.hdf5"
fid  = h5open(file)
runs = keys(fid)
run  = first(runs)


N_replicas = read(fid[run],"N_replicas")
N_repeats = read(fid[run],"N_repeats")
a_last = zeros(N_replicas,N_repeats)
p_last = zeros(N_replicas,N_repeats)
for i in 0:N_replicas-1, j in 0:N_repeats-1
    a_last[i+1,j+1] = fid[run]["$j/Rep_$i/a"][end]
    p_last[i+1,j+1] = fid[run]["$j/Rep_$i/plaq"][end]
end

a0  = dropdims(mean(a_last,dims=2),dims=2)
Δa0 = dropdims(std(a_last,dims=2),dims=2)/sqrt(N_repeats)
p0  = dropdims(mean(p_last,dims=2),dims=2)

scatter(p0,a0,yerr=Δa0)