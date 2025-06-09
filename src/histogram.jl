# In python, we are using the following precision
# Mpmath settings:
#   mp.prec = 53                [default: 53]
#   mp.dps = 15                 [default: 15]
#   mp.trap_complex = False     [default: False]
setprecision(BigFloat, 53)

function log_partition_function(a, S, beta)
    # David uses a different sign for a
    dS     = S[2] - S[1]    
    pi_exp = 0.0
    Z      = BigFloat(0)
    for (ai, Si) in zip(a,S)
        A = beta - ai
        exp_factor  = exp(BigFloat(pi_exp + Si*beta - ai*dS/2))
        sinh_factor = iszero(A) ? dS/2 : sinh(A*dS/2)/A
        Z      += exp_factor*sinh_factor
        pi_exp -= ai*dS
    end
    return Float64(log(2Z))
end
_E_in_interval(E,S0,dS) = E >= (S0 - dS/2) && E < (S0 + dS/2)
function log_rho(E, S, dS, a; cumsum_a = cumsum(a))
    for i in eachindex(S, a)   
        if _E_in_interval(E,S[i],dS)
            log_ρ = a[i] * (S[i] - dS/2 - E) - cumsum_a[i-1]*dS
            return log_ρ
        end
    end
end
function _set_up_histogram(fid,run)
    a, S_all = LLRParsing.a_vs_central_action_repeats(fid,run;ind=nothing)[1:2]
    Nl   = read(fid[run],"Nl")
    Nt   = read(fid[run],"Nt")
    S    = unique(S_all) 
    V    = Nt*Nl^3
    return a, S, Nt, Nl, V
end
function probability_density_repeats(fid, run, beta; kws...)
    a, S, Nt, Nl, V = _set_up_histogram(fid,run)
    return probability_density_repeats(a, S, beta, V; kws...)
end
function probability_density_repeats(a, S, beta, V; nbins=length(S))
    dS  = S[2] - S[1]
    E   = range(minimum(S) + dS, maximum(S), length=nbins)
    ups = @. E/(6V) 
    repeats = last(size(a)[2])
    P = zeros(nbins,repeats)

    for i in 1:repeats
        logZ = log_partition_function(a[:,i], S, beta)
        csa  = cumsum(a[:,i])
        for j in 1:nbins
            log_ρ  = log_rho(E[j], S, dS, a[:,i]; cumsum_a=csa)
            P[j,i] = exp(log_ρ + beta*E[j] - logZ)
        end
    end
    return ups, P, V, dS
end
function probability_density(fid, run, beta; kws...)
    ups, prob, V, dS = probability_density_repeats(fid, run, beta; kws...)
    P  = dropdims(mean(prob,dims=2),dims=2)
    ΔP = dropdims(std(prob,dims=2),dims=2)/sqrt(size(prob)[2])
    covP = cov(prob,dims=2)/size(prob)[2] 
    return ups, P, ΔP, covP, V, dS
end
function probability_density(a, S, beta, V; kws...)
    ups, prob, V, dS = probability_density_repeats(a, S, beta, V; kws...)
    P  = dropdims(mean(prob,dims=2),dims=2)
    ΔP = dropdims(std(prob,dims=2),dims=2)/sqrt(size(prob)[2])
    covP = cov(prob,dims=2)/size(prob)[2]
    return ups, P, ΔP, covP, V, dS
end
function plot_plaquette_histogram!(plt,fid,run,beta;kws...)
    a, S, Nt, Nl, V = _set_up_histogram(fid,run)
    ups, P, ΔP, covP, V, dS = probability_density(a, S, beta, V)
    label  = "$(Nt)x$(Nl): ΔE=$(round(2(dS)/6V,sigdigits=1))"
    xlabel = L"u_p"
    ylabel = L"P_{\beta}(u_p)"
    plot!(plt;xlabel,ylabel,yticks=:none,left_margin=5Plots.mm)
    plot!(plt,ups,P*6V;label,ribbon=ΔP*6V,kws...)
end