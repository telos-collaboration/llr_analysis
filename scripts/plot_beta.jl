using DelimitedFiles
using Plots
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=10,labelfontsize=12,left_margin=1Plots.mm)

function read_critical_betas(file)
    data = readdlm(file,',',skipstart=1)
    βc   = data[:,4]
    Δβc  = data[:,5]
    N = size(data,1)
    βc11  = βc[1:N÷2]
    Δβc11 = Δβc[1:N÷2]
    βc22  = βc[N÷2+1:end]
    Δβc22 = Δβc[N÷2+1:end]
    return βc11, Δβc11, βc22, Δβc22
end

L = [48,48,56,56,56,64,72,80]
file = "data_assets/critical_beta_Nt5.csv"
βc11, Δβc11, βc22, Δβc22 = read_critical_betas(file)

plt1 = plot(;title=L"N_t = 5",ylabel=L"\beta_c",xlabel=L"spatial volume $N_l$")
scatter!(plt1,L,βc11,yerr=Δβc11,label=L"$P_\beta$ peak ratio 1:1")
scatter!(plt1,L,βc22,yerr=Δβc22,label=L"$P_\beta$ peak ratio 2:1")

L = [20,24,28,40,48]
file = "data_assets/critical_beta_Nt4.csv"
βc11, Δβc11, βc22, Δβc22 = read_critical_betas(file)

plt2 = plot(;title=L"N_t = 4",ylabel=L"\beta_c",xlabel=L"spatial volume $N_l$")
scatter!(plt2,L,βc11,yerr=Δβc11,label=L"$P_\beta$ peak ratio 1:1")
scatter!(plt2,L,βc22,yerr=Δβc22,label=L"$P_\beta$ peak ratio 2:1")

savefig(plt1,"assets/plots/critical_beta_volumes_Nt5.pdf")
savefig(plt2,"assets/plots/critical_beta_volumes_Nt4.pdf")