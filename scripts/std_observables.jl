using Pkg; Pkg.activate(".")
using BenchmarkTools
using HDF5
using StatsBase

# For higher moments I use the didcated functions because they are faster
# I have checked that the numerical precision is sufficient
binder_cumulant(x) = 1 - moment(x,4,0.0) / mean(f->f^2,x)^2 / 3
poly_susceptibility(x,Nl) = var(abs.(x))*Nl^3
specific_heat_plaq(x,Nl,Nt) = var(x)*6Nl^3*Nt

# TODO: Give variance of the mean of those statistical uncertainties analytically
function std_observables(f,ens)
    plaq = f[ens]["plaquette"][]
    poly = f[ens]["polyakov_loop"][]
    Nt = f[ens]["Nt"][]
    Nl = f[ens]["Nl"][]

    binder_plaq = binder_cumulant(plaq)
    poly_sus = poly_susceptibility(poly,Nl)
    sh_plaq  = specific_heat_plaq(plaq,Nl,Nt)
    poly_vev = mean(abs,poly)
    plaq_vev = mean(plaq)
    return binder_plaq, poly_sus, sh_plaq, poly_vev, plaq_vev, plaq, poly
end

h5file = "test.hdf5"
f      = h5open(h5file)
ens    = "ImportanceSampling/4x20/7.32/"

binder_plaq, poly_sus, sh_plaq, poly_vev, plaq_vev, plaq, poly = std_observables(f,ens)

using MadrasSokal
madras_sokal_time(plaq)
madras_sokal_time(abs.(poly))

# generate a resample of the original correlator matrix
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

    return resampled
end
