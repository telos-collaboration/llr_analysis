using Pkg; Pkg.activate(".")
using HDF5
using LLRParsing
using Plots
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
betas = [7.48967, 7.48970, 7.48982, 7.48969, 7.48975]
runs  = keys(fid)

function histogram_comparison(fid,runs,betas)
    plt = plot(xlabel=L"u_p",ylabel=L"P_{\beta}(u_p)",yticks=:none,left_margin=5Plots.mm)
    for (i,run) in enumerate(runs)
        @show run
        Nl   = read(fid[run],"Nl")
        Nt   = read(fid[run],"Nt")
        xl   = (0.5885,0.5905)
        beta = betas[i]

        ups,P,ΔP,V,dS = probability_density(fid, run, beta)
        label = "$(Nt)x$(Nl): ΔE=$(round(2(dS)/6V,sigdigits=1))"
        plot!(plt,ups,P*6V;label,ribbon=ΔP*6V,xlims=xl)
    end
    return plt
end
plt = histogram_comparison(fid,runs,betas)
display(plt)