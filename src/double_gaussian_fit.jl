Gaussian(x,A,μ,σ) = A * exp(-(x-μ)^2 / σ^2 / 2)
DoubleGaussian(x,A1,μ1,σ1,A2,μ2,σ2) = Gaussian(x,A1,μ1,σ1) + Gaussian(x,A2,μ2,σ2)
@. modelDG(x, p) = DoubleGaussian(x,p[1],p[2],p[3],p[4],p[5],p[6])
function initial_param_double_gaussian(ups,P;w=5)
    δups = ups[2]-ups[1]
    pks  = findmaxima(P,w)
    pks  = peakproms(pks)
    pks  = peakwidths(pks)
    p0   = [ pks.heights[1],ups[pks.indices[1]],pks.widths[1]*δups,
    pks.heights[2],ups[pks.indices[2]],pks.widths[2]*δups ]
    return p0
end
function fitting_range_double_gaussian(P;w=5,N=0.5)
    # determine fitting range: 
    # Include half the points between peak and minimum 
    pks = findmaxima(P,w)
    mns = findminima(P,w)
    δ1  = abs(pks.indices[1]-mns.indices[1])*N
    δ2  = abs(pks.indices[2]-mns.indices[1])*N
    δ   = Int(round((δ1+δ2)/2))
    range1 = pks.indices[1]-δ:pks.indices[1]+δ
    range2 = pks.indices[2]-δ:pks.indices[2]+δ
    range  = vcat(collect(range1),collect(range2)) 
    return range
end
function fit_double_gaussian(ups::Vector{T},P::Vector{T},covP::Matrix{T};kws...) where T
    p0  = initial_param_double_gaussian(ups,P;kws...)
    r   = fitting_range_double_gaussian(P; kws...)
    fit = curve_fit(modelDG, ups[r], P[r], p0)
    return fit 
end
function fit_double_gaussian(ups::Vector{T},P::Vector{T};kws...) where T
    p0  = initial_param_double_gaussian(ups,P;kws...)
    r   = fitting_range_double_gaussian(P; kws...)
    fit = curve_fit(modelDG, ups[r], P[r], p0)
    return fit 
end
function fit_double_gaussian(fid::HDF5.File, run, beta;kws...)
    ups,P,ΔP,covP,V,dS = probability_density(fid, run, beta)
    fit = fit_double_gaussian(ups,P,covP;kws...)
    return fit
end
function plot_double_gaussian_fit!(plt,fid, run, beta;kws...)
    ups,P,ΔP,covP,V,dS = probability_density(fid, run, beta)
    fit = fit_double_gaussian(ups,P,covP;kws...)
    fitted = 6 .* V .* modelDG(ups,fit.param)
    plot!(plt,ups,fitted,label="double Gaussian fit",lw=3)
    return plt
end
