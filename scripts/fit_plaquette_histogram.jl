using HDF5
using LLRParsing
using Plots
using Statistics
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
betas = [7.48967, 7.48970, 7.48982, 7.48969, 7.48975]
runs  = keys(fid)
ind   = 5
beta  = betas[ind]
run   = runs[ind] 

xl   = (0.5885,0.5905)
ups,P,ΔP,V,dS = probability_density(fid, run, beta)
plot_plaquette_histogram!(plt,fid,run,betas[i];xlims=xl)

using Peaks
δups = ups[2]-ups[1]
pks  = findmaxima(vec(P),50)
pks  = peakproms(pks)
pks  = peakwidths(pks)

using LsqFit
Gaussian(x,A,μ,σ) = A * exp(-(x-μ)^2 / σ^2 / 2)
DoubleGaussian(x,A1,μ1,σ1,A2,μ2,σ2) = Gaussian(x,A1,μ1,σ1) + Gaussian(x,A2,μ2,σ2)
@. model(x, p) = DoubleGaussian(x,p[1],p[2],p[3],p[4],p[5],p[6])
p0  = [pks.heights[1],ups[pks.indices[1]],pks.widths[1]*δups,pks.heights[2],ups[pks.indices[2]],pks.widths[2]*δups]
fit = curve_fit(model, ups, vec(P), p0)
y0  = model(ups,p0)
y   = model(ups,fit.param)
plot!(plt,ups,y,label="double Gaussian fit",lw=3)