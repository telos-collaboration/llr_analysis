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
function fit_double_gaussian(fid, run, beta)
    ups,P,ΔP,V,dS = probability_density(fid, run, beta)
    fit = fit_double_gaussian(ups,P)
    return fit
end
function plot_double_gaussian_fit!(plt,fid, run, beta)
    ups,P,ΔP,V,dS = probability_density(fid, run, beta)
    fit = fit_double_gaussian(ups,P)
    fitted = 6 .* V .* modelDG(ups,fit.param)
    plot!(plt,ups,fitted,label="double Gaussian fit",lw=3)
    return plt
end

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
betas = [7.48967, 7.48970, 7.48982, 7.48969, 7.48975]
runs  = keys(fid)

for i in eachindex(betas,runs)
    beta = betas[i]
    run  = runs[i] 

    plt = plot()
    plot_plaquette_histogram!(plt,fid,run,beta;xlims=(0.5885,0.5905))
    plot_double_gaussian_fit!(plt,fid,run,beta)
    display(plt)
end