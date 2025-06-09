using LLRParsing
using HDF5
using Statistics
using Plots

function jackknife_resamples(obs)
    nops,nconf = size(obs)
    samples = similar(obs)
    # Strategy: Sum over all configs then remove one
    tmp = dropdims(sum(obs,dims=2),dims=2)
    for index in 1:nconf    
        @. samples[:,index] = (tmp - obs[:,index])/(nconf-1) 
    end
    return samples
end
function histogram_jackknife_fit(fid,run)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,run)
    a_jk   = jackknife_resamples(a) 
    fitted = similar(a_jk)
    ups    = similar(S)
    for i in axes(a_jk,2)
        ai   = a_jk[:,i:i]
        beta = LLRParsing.beta_at_equal_heights(ai, S, V)
        ups, P, ΔP, covP, V, dS = LLRParsing.probability_density(ai, S, beta, V)
        fit  = LLRParsing.fit_double_gaussian(ups,P)
        data = 6 .* V .* LLRParsing.modelDG(ups,fit.param)
        fitted[:,i] = data 
    end
    N  = size(fitted)[2]
    f  = dropdims(mean(fitted,dims=2),dims=2)
    Δf = sqrt(N-1)*dropdims(std(fitted,dims=2),dims=2)
    return ups, f, Δf
end

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
runs  = keys(fid)
run   = runs[end-2]

ups, f, Δf = histogram_jackknife_fit(fid,run)
plot(ups,f,ribbon=Δf)