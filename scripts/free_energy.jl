using LLRParsing
using HDF5
using Statistics
using Plots
using LaTeXStrings
using LinearAlgebra
using PCHIPInterpolation
using Peaks
using Roots
using ArgParse
gr(fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=1Plots.mm)

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
            s[i,r] = (logρ + logρ0)/V
            t[i,r] = 1/a0[i]
            f[i,r] = u - t[i,r]*s[i,r]
        end
    end
    return t, f, s
end
function thermodynamic_potentials(a, S0, V; kws...)
    t_r, f_r, s_r = thermodynamic_potentials_repeats(a, S0, V; kws...)
    N  = size(t_r,2)
    t  = dropdims(mean(t_r,dims=2),dims=2)
    f  = dropdims(mean(f_r,dims=2),dims=2)
    s  = dropdims(mean(s_r,dims=2),dims=2)
    Δt = dropdims(std(t_r,dims=2),dims=2) ./ sqrt(N)
    Δf = dropdims(std(f_r,dims=2),dims=2) ./ sqrt(N)
    Δs = dropdims(std(s_r,dims=2),dims=2) ./ sqrt(N)
    return t, Δt, f, Δf, s, Δs
end
function plot_free_energies(file,plotdir)
    h5dset = h5open(file)
    runs   = keys(h5dset)
    tmin, tmax = +Inf, -Inf
    fmin, fmax = -Inf, +Inf
    plt_s = plot()
    for r in runs
        a, Δa, S0, _ = a_vs_central_action(h5dset,r)
        t, Δt, f, Δf, s, Δs = thermodynamic_potentials(h5dset,r)
        pks = only(findmaxima(a,5).indices)
        mns = only(findminima(a,5).indices)
        r1  = 1:pks
        r2  = mns:length(a)   

        perm1 = sortperm(t[r1])
        perm2 = sortperm(t[r2])
        itp1 = Interpolator(t[r1][perm1],f[r1][perm1])
        itp2 = Interpolator(t[r2][perm2],f[r2][perm2])

        g(x)  = itp1(x) - itp2(x)
        t1,t2 = extrema(filter(t -> isfinite(g(t)), vcat(t[r1],t[r2]) ))
        tc    = find_zero(g,(t1,t2))
        f     = (f .- itp1(tc)) 
        @. f  = f * 10^6
        @. Δf = Δf * 10^6
        
        # set plot limits
        tmin = min(t[pks],tmin) 
        tmax = max(t[mns],tmax) 
        fmin = max(f[pks],fmin) 
        fmax = min(f[mns],fmax)
        δt, δf = tmax-tmin, fmax-fmin 
        
        plt = plot(title=LLRParsing.fancy_title(r)*" - no entropy subtraction")
        plot!(;ylabel=L"(f - f_c^+ )/ 10^{-6}", xlabel=L"t = 1/a_n")
        plot!(plt,t,f,xerr=Δt,yerr=Δf,ms=1,label="")
        plot!(plt,ylims=(fmin-δf,fmax+δf/2),xlims=(tmin-δt/4,tmax+δt/4))
        ispath(plotdir) || mkpath(plotdir)
        savefig(plt,joinpath(plotdir,"$r.pdf"))
        plot!(plt_s,t,s,xerr=Δt,yerr=Δs,ms=1,label=r)
    end
    Nt = read(h5dset[first(runs)],"Nt")
    plot!(legend=:bottomright, xlabel=L"t = 1/a_n", ylabel=L"unsubtracted entropy $s = \log(\rho)$")
    savefig(plt_s,joinpath(plotdir,"entropy_Nt$Nt.pdf"))
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
            help = "HDF5 file containing the sorted results"
            required = true
        "--plot_dir"
            help = "Where to save the plots"
            required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file    = args["h5file"]
    plotdir = args["plot_dir"]
    plot_free_energies(file, plotdir)
end
main()