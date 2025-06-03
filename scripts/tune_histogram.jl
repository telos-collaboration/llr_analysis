using HDF5
using LLRParsing
using Plots
using Statistics
using LaTeXStrings
using Peaks
using Roots
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function peak_height_difference(fid,run,β;w=50)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,run)
    return peak_height_difference(a, S, β, V;w)
end
function peak_height_difference(a, S, β, V;w=50)
    ups, P, ΔP, V, dS = probability_density(a, S, β, V; nbins=1000)
    pks = findmaxima(P,w)
    n_peaks = length(pks.indices)
    if n_peaks == 2
        heightdiff = pks.heights[2] - pks.heights[1]
        return heightdiff
    else 
        @warn "Found $n_peaks peak(s)"
        return nothing
    end
end
function bracket_critical_beta(fid,run;w=50,Nint=50)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,run)
    βl, βu = extrema(a)
    βold = (βl + βu)/2
    βint = abs(βl - βu)/2
    for _ in 1:Nint
        diff    = peak_height_difference(a, S, βold, V; w)
        βnew    = diff>0 ? βold - βint/Nint : βold + βint/Nint 
        newdiff = peak_height_difference(a, S, βnew, V; w)
        if sign(diff) != sign(newdiff)
            # use extrame function to sort bracket of equal height
            return extrema((βold,βnew))
        end
        # move to the next interval
        βold, diff = βnew, newdiff
    end
end
function beta_at_equal_heights(fid,run)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,run)
    βl, βu = bracket_critical_beta(fid,run)
    f(β)   = peak_height_difference(a, S, β, V)
    βc     = find_zero(f, (βl, βu), Bisection())
    return βc, βl, βu
end

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
runs  = keys(fid)
for run in runs
    @show runs
    βc, βl, βu = beta_at_equal_heights(fid,run)
    @show βc, βl, βu
end 
@time beta_at_equal_heights(fid,last(runs))
@profview beta_at_equal_heights(fid,last(runs))
βc, βl, βu = beta_at_equal_heights(fid,last(runs))