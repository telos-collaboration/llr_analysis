using Pkg; Pkg.activate(".")
using BenchmarkTools
using HDF5
using StatsBase
using Statistics
using MadrasSokal

# For higher moments I use the didcated functions because they are faster
# I have checked that the numerical precision is sufficient
binder_cumulant(x) = 1 - moment(x,4,0.0) / mean(f->f^2,x)^2 / 3
poly_susceptibility(x,Nl) = var(abs.(x))*Nl^3
specific_heat_plaq(x,Nl,Nt) = var(x)*6Nl^3*Nt

# TODO: Give variance of the mean of those statistical uncertainties analytically
#       Bin data instead of separating the measurements by τint
function jackknife_resample_1d_reduction(x,f)
    resampled = similar(x)
    tmp = zeros(eltype(x),(length(x)-1))
    for index in eachindex(x)    
        for i in eachindex(x)
            index == i && continue
            j = i < index ? i : i - 1
            tmp[j] = x[i]
        end
        resampled[index] = f(tmp)
    end
    return apply_jackknife(resampled)
end
function apply_jackknife(obs::AbstractVector)
    N  = length(obs)
    O  = mean(obs)
    ΔO = sqrt(N-1)*std(obs,corrected=false)
    return O, ΔO
end
function std_observables(f,ens)
    plaq = f[ens]["plaquette"][]
    poly = f[ens]["polyakov_loop"][]
    Nt = f[ens]["Nt"][]
    Nl = f[ens]["Nl"][]

    binder_plaq = binder_cumulant(plaq)
    sh_plaq  = specific_heat_plaq(plaq,Nl,Nt)
    plaq_vev = mean(plaq)
    poly_sus = poly_susceptibility(poly,Nl)
    poly_vev = mean(abs,poly)

    # Estimate autocorrelation times: 
    # To speed it up I first obtain the exponential autocorrelation time 'τexp' and then
    # compute the integrated autocorrelation time 'τint' using the Madras-Sokal windowing
    # technique on measurements seperated by 'τexp'.
    τexp_plaq = exponential_autocorrelation_time(plaq)
    step_plaq = Int(ceil(τexp_plaq))
    τint_plaq, Δτint_plaq = step_plaq.*madras_sokal_time(plaq[1:step_plaq:end])

    τexp_poly = exponential_autocorrelation_time(abs.(poly))
    step_poly = Int(ceil(τexp_poly/4))
    τint_poly, Δτint_poly = step_poly.*madras_sokal_time(abs.(poly[1:step_poly:end]))

    # Determine uncertainties from data seperated by τint
    plaq_τ =      plaq[1:Int(ceil(τint_plaq)):end]
    poly_τ = abs.(poly[1:Int(ceil(τint_poly)):end])

    Δbinder_plaq = jackknife_resample_1d_reduction(plaq_τ,x -> binder_cumulant(x))[2]
    Δsh_plaq  = jackknife_resample_1d_reduction(plaq_τ,x -> specific_heat_plaq(x,Nl,Nt))[2]
    Δplaq_vev = std(plaq_τ)/sqrt(length(plaq_τ))
    Δpoly_sus = jackknife_resample_1d_reduction(poly_τ,x -> poly_susceptibility(x,Nl))[2]
    Δpoly_vev = std(poly_τ)/sqrt(length(poly_τ))
    
    return binder_plaq, Δbinder_plaq, sh_plaq, Δsh_plaq, plaq_vev, Δplaq_vev, poly_sus, Δpoly_sus, poly_vev, Δpoly_vev
end

h5file = "LLR_data.hdf5"
f      = h5open(h5file)
ens    = "ImportanceSampling/4x20/7.32/"
binder, Δbinder, sh_plaq, Δsh_plaq, plaq, Δplaq, poly_sus, Δpoly_sus, poly, Δpoly = std_observables(f,ens)