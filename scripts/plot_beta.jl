using DelimitedFiles
using Plots

file = "data_assets/critical_beta_Nt5.csv"
data = readdlm(file,',',skipstart=1)
βc   = data[:,4]
Δβc  = data[:,5]

N = size(data,1)

βc11  = βc[1:N÷2]
Δβc11 = Δβc[1:N÷2]
βc22  = βc[N÷2+1:end]
Δβc22 = Δβc[N÷2+1:end]

L = [20,24,28,40,48]
L = [48,48,56,56,56,64,72,80]

scatter(L,βc11,yerr=Δβc11)
scatter!(L,βc22,yerr=Δβc22)