using HDF5
using LLRParsing
using Plots
using Statistics
using LaTeXStrings
using Peaks
using LsqFit
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

Gaussian(x,A,μ,σ) = A * exp(-(x-μ)^2 / σ^2 / 2)
DoubleGaussian(x,A1,μ1,σ1,A2,μ2,σ2) = Gaussian(x,A1,μ1,σ1) + Gaussian(x,A2,μ2,σ2)
@. modelDG(x, p) = DoubleGaussian(x,p[1],p[2],p[3],p[4],p[5],p[6])
function initial_param_double_gaussian(ups,P)
    δups = ups[2]-ups[1]
    pks  = findmaxima(vec(P),50)
    pks  = peakproms(pks)
    pks  = peakwidths(pks)
    p0   = [ pks.heights[1],ups[pks.indices[1]],pks.widths[1]*δups,
    pks.heights[2],ups[pks.indices[2]],pks.widths[2]*δups ]
    return p0
end
function fit_double_gaussian(ups,P)
    p0  = initial_param_double_gaussian(ups,P)
    fit = curve_fit(modelDG, ups, vec(P), p0)
    return fit 
end

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
betas = [7.48967, 7.48970, 7.48982, 7.48969, 7.48975]
runs  = keys(fid)
ind   = 5
beta  = betas[ind]
run   = runs[ind] 

ups,P,ΔP,V,dS = probability_density(fid, run, beta)
fit = fit_double_gaussian(ups,P)
plt = plot_plaquette_histogram!(plot(),fid,run,betas[ind];xlims=(0.5885,0.5905))
y   = modelDG(ups,fit.param)
plot!(plt,ups,y*6V,label="double Gaussian fit",lw=3)