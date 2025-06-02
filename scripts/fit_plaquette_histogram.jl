using HDF5
using LLRParsing
using Plots
using Statistics
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
beta  = 7.48975
run   = keys(fid)[end]

Nl   = read(fid[run],"Nl")
Nt   = read(fid[run],"Nt")
xl   = (0.5885,0.5905)

ups,P,ΔP,V,dS = probability_density(fid, run, beta)

plt   = plot(xlabel=L"u_p",ylabel=L"P_{\beta}(u_p)",yticks=:none,left_margin=5Plots.mm)
label = "$(Nt)x$(Nl): ΔE=$(round(2(dS)/6V,sigdigits=1))"
plot!(plt,ups,P;label,ribbon=ΔP,xlims=xl,marker=:circle,ms=2,markeralpha=0.8)


using Peaks
δups = ups[2]-ups[1]
pks  = findmaxima(vec(P),50)
pks  = peakproms(pks)
pks  = peakwidths(pks)

pks.indices

using LsqFit
Gaussian(x,A,μ,σ) = A * exp(-(x-μ)^2 / σ^2 / 2)
DoubleGaussian(x,A1,μ1,σ1,A2,μ2,σ2) = Gaussian(x,A1,μ1,σ1) + Gaussian(x,A2,μ2,σ2)
@. model(x, p) = DoubleGaussian(x,p[1],p[2],p[3],p[4],p[5],p[6])
p0  = [pks.heights[1],ups[pks.indices[1]],pks.widths[1]*δups,pks.heights[2],ups[pks.indices[2]],pks.widths[2]*δups]
fit = curve_fit(model, ups, vec(P), p0)
y0  = model(ups,p0)
y   = model(ups,fit.param)
plot!(plt,ups,y,label="double Gaussian fit",lw=3)