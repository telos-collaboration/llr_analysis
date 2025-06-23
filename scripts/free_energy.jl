using LLRParsing
using HDF5
using Statistics
using Plots
using LaTeXStrings
using LinearAlgebra
using PCHIPInterpolation
using Peaks
using Roots
gr(fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function thermodynamic_potentials_repeats(h5dset,run; kws...)
    a, S0_all, ind = LLRParsing.a_vs_central_action_repeats(h5dset,run)
    Nl = read(h5dset[run],"Nl")
    Nt = read(h5dset[run],"Nt")
    V  = Nt*Nl^3
    S0 = dropdims(unique(S0_all,dims=2),dims=2)
    return thermodynamic_potentials_repeats(a, S0, V; kws...)
end
function thermodynamic_potentials(h5dset,run; kws...)
    a, S0_all, ind = LLRParsing.a_vs_central_action_repeats(h5dset,run)
    Nl = read(h5dset[run],"Nl")
    Nt = read(h5dset[run],"Nt")
    V  = Nt*Nl^3
    S0 = dropdims(unique(S0_all,dims=2),dims=2)
    return thermodynamic_potentials(a, S0, V; kws...)
end
function thermodynamic_potentials_repeats(a, S0, V; logρ0 = 0.0)
    dS = S0[2] - S0[1]
    E  = range(minimum(S0) + dS, maximum(S0), length=length(S0))
    s, t, f, u = zero(a),zero(a),zero(a),zero(a)
    replicas = size(a,1)
    repeats  = size(a,2)
    for r in 1:repeats
        a0 = a[:,r]
        cs = cumsum(a0)
        for i in 1:replicas
            logρ  = LLRParsing.log_rho(E[i], S0, dS, a0; cumsum_a=cs)
            u      = (6V - E[i])/V
            s      = (logρ + logρ0)/V
            t[i,r] = 1/a0[i]
            f[i,r] = u - t[i,r]*s
        end
    end
    return t, f
end
function thermodynamic_potentials(a, S0, V; kws...)
    t_r, f_r = thermodynamic_potentials_repeats(a, S0, V; kws...)
    N  = size(t_r,2)
    t  = dropdims(mean(t_r,dims=2),dims=2)
    f  = dropdims(mean(f_r,dims=2),dims=2)
    Δt = dropdims(std(t_r,dims=2),dims=2) ./ sqrt(N)
    Δf = dropdims(std(f_r,dims=2),dims=2) ./ sqrt(N)
    return t, Δt, f, Δf
end
function main()
    file   = "data_assets/Sp4_Nt5_sorted.hdf5"
    h5dset = h5open(file)
    runs   = keys(h5dset)
    tmin, tmax = +Inf, -Inf
    fmin, fmax = -Inf, +Inf
    for r in runs
        plt    = plot()
        a, Δa, S0, _ = a_vs_central_action(h5dset,r)
        t, Δt, f, Δf = thermodynamic_potentials(h5dset,r)
        pks = only(findmaxima(a,5).indices)
        mns = only(findminima(a,5).indices)
        r1  = 1:pks
        r2  = mns:length(a)   

        perm1 = sortperm(t[r1])
        perm2 = sortperm(t[r2])
        itp1 = Interpolator(t[r1][perm1],f[r1][perm1])
        itp2 = Interpolator(t[r2][perm2],f[r2][perm2])

        g(x)   = itp1(x) - itp2(x)
        t1, t2 = extrema(filter(t -> isfinite(g(t)), vcat(t[r1],t[r2]) ))
        tc     = find_zero(g,(t1,t2))
        f      = f .- itp1(tc)
        plot!(plt,t,f,xerr=Δt,yerr=Δf,label="$r",ms=1)

        # set plot limits
        tmin = min(t[pks],tmin) 
        tmax = max(t[mns],tmax) 
        fmin = max(f[pks],fmin) 
        fmax = min(f[mns],fmax)
        δt, δf = tmax-tmin, fmax-fmin 
        plot!(plt,ylims=(fmin-δf,fmax+δf),xlims=(tmin-δt/4,tmax+δt/4))
        display(plt)
    end
    return plt
end
plt = main()