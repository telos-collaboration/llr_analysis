# In python, we are using the following precision
# Mpmath settings:
#   mp.prec = 53                [default: 53]
#   mp.dps = 15                 [default: 15]
#   mp.trap_complex = False     [default: False]

function log_partition_function(a, S, beta)
    
    # David uses a different sign for a
    a      = -a
    dS     = S[2] - S[1]    
    pi_exp = BigFloat(0)
    Z      = BigFloat(0)
    
    for (ai, Si) in zip(a,S)
        A = ai + beta
        if iszero(A)
            Z += exp(pi_exp - Si*ai + ai*dS/2)  * dS/2
        else
            full_exp    = exp(pi_exp + Si*beta + ai*dS/2)
            sinf_factor = sinh(A*dS/2)
            Z += full_exp*sinf_factor/A
        end
        pi_exp += ai*dS
    end
    return Float64(log(2Z))
end
function log_rho(E, S, dS, a)
    # David uses a different sign for a
    a = -a
    S_shifted = @. S - dS/2
    pi_exp = BigFloat(0)
    log_ρ  = BigFloat(0)
    for (Si,ai) in zip(S_shifted, a)   
        if E >= Si && E < (Si + dS)
            log_ρ = (ai * (E - Si)) + pi_exp
            break
        else
            pi_exp += ai*dS
        end
    end
    @assert !iszero(log_ρ)
    return Float64(log_ρ)
end
function probability_density(a, S, beta, V; nbins=1000)
    up   = S/(6V)
    dS   = S[2] - S[1]
    δup  = dS/(6V)
    ups  = range(minimum(up) + δup, maximum(up), length=nbins)
    E    = @. ups*6V 

    repeats = last(size(a))
    probability_density = zeros(nbins,repeats)
    
    for i in 1:repeats
        logZ = log_partition_function(a[:,i], S, beta)
        log_ρ = log_rho.(E, Ref(S), dS, Ref(a[:,i]))
        probability_density[:,i] = @. exp(log_ρ + beta*ups*V*6 - logZ)
    end
    return ups, probability_density
end
function probability_density_repeats(fid, run, beta; kws...)
    a, S_all = a_vs_central_action_repeats(fid,run;ind=nothing)[1:2]
    Nl   = read(fid[run],"Nl")
    Nt   = read(fid[run],"Nt")
    V    = Nl^3 * Nt
    S    = unique(S_all) 
    dS   = S[2] - S[1]
    
    ups,P = probability_density(a, S, beta, V; kws...)
    return ups, P, V, dS
end
function probability_density(fid, run, beta; kws...)
    ups, prob, V, dS = probability_density_repeats(fid, run, beta; kws...)
    P  = dropdims(mean(prob,dims=2),dims=2)
    ΔP = dropdims(std(prob,dims=2),dims=2)/sqrt(size(prob)[2])
    return ups, P, ΔP, V, dS
end
function plot_plaquette_histogram!(plt,fid,run,beta;kws...)
    Nl  = read(fid[run],"Nl")
    Nt  = read(fid[run],"Nt")
    ups,P,ΔP,V,dS = probability_density(fid, run, beta)
    label  = "$(Nt)x$(Nl): ΔE=$(round(2(dS)/6V,sigdigits=1))"
    xlabel = L"u_p"
    ylabel = L"P_{\beta}(u_p)"
    plot!(plt;xlabel,ylabel,yticks=:none,left_margin=5Plots.mm)
    plot!(plt,ups,P*6V;label,ribbon=ΔP*6V,kws...)
end