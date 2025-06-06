using HDF5
using LLRParsing
using Plots
using Statistics
using LaTeXStrings
using Peaks
using Roots
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)
"""
    peak_height_difference(args;w=50)

    Calculate the height difference between the two peaks in the histogram.

    `w` is the minimal number of intervals that the peaks need to be separated
    to count as proper peaks. This variable probably needs to be tuned depending
    on the noise in the dataset.
"""
function peak_height_difference(fid,run,β;w=50)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,run)
    return peak_height_difference(a, S, β, V;w)
end
function peak_height_difference(a, S, β, V;w=50)
    ups, P, ΔP, V, dS = probability_density(a, S, β, V)
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
"""
    bracket_critical_beta(args; w=50, Nint=50)

    Find two value of beta for which the height difference has a different sign.
    The two values can then be used to determine the critical beta using a Bisection
    search to find the value for which the height difference is vanishing.

    `w` is the minimal required separation of the peaks.

    Here I first look for an initial β for which we have to peaks. Then, depending on 
    the sign of the height difference, β is lowered or raised in fixed steps. The maximal 
    number of steps is given by `Nint`.    
"""
function bracket_critical_beta(fid,run; kws...)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,run)
    bracket_critical_beta(a, S, V ; kws...)
end
function bracket_critical_beta(a, S, V; kws...)
    βmin, βmax = extrema(a)
    β0     = (βl + βu)/2
    bracket_critical_beta(a, S, V, β0, βmin, βmax; kws...)
end
function bracket_critical_beta(a, S, V, β0, βmin, βmax;w=50,Nint=50)
    βold = β0
    βint = abs(βmax - βmin)/2
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
"""
    Determine the critical coupling by first finding an interval for which the height
    difference of the histogram peaks vanish
"""
function beta_at_equal_heights(fid,run;kws...)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,run)
    return beta_at_equal_heights(a, S, V;kws...)
end
function beta_at_equal_heights(a, S, V; kws...)
    βl, βu = bracket_critical_beta(a, S, V;kws...)
    f(β)   = peak_height_difference(a, S, β, V)
    βc     = find_zero(f, (βl, βu), Bisection())
    return βc, βl, βu
end
function beta_at_equal_heights(a, S, V, β0, βmin, βmax; kws...)
    βl, βu = bracket_critical_beta(a, S, V, β0, βmin, βmax;kws...)
    f(β)   = peak_height_difference(a, S, β, V)
    βc     = find_zero(f, (βl, βu), Bisection())
    return βc, βl, βu
end
function count_peaks(a, S, β, V; w = 50, kws...)
    ups, P, ΔP, V, dS = probability_density(a, S, β, V)
    pks = findmaxima(P,w)
    n_peaks = length(pks.indices)
    return n_peaks
end
function find_initial_double_peak(a,S,β0,βmin,βmax,V;Nint=100,w=50)
    # start looking at β0 in the interval [βmin,βmax]
    np = count_peaks(a, S, β0, V; w)
    if np == 2
        return β0
    end
    # If there is no double peak at β0, then partition the interval
    # into Nint intervals and check them in order of proximity to βmin.
    βs          = collect(range(βmin,βmax,length=Nint))
    ordered_ind = sort(eachindex(βs), by = ind -> abs(βs[ind] - β0))
    for i in ordered_ind
        np = count_peaks(a, S, βs[i], V; w)
        if np == 2
            return βs[i]
        end
    end
end

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
runs  = keys(fid)
#for run in runs
#    @show runs
#    βc, βl, βu = beta_at_equal_heights(fid,run)
#    @show βc, βl, βu
#end 

#@time beta_at_equal_heights(fid,last(runs))
#@profview beta_at_equal_heights(fid,last(runs))

for run in runs
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,run)
    βmin, βmax      = extrema(a)
    β0              = (βmin + βmax)/2
    βdg             = find_initial_double_peak(a,S,β0,βmin,βmax,V)
    βl, βu          = bracket_critical_beta(a, S, V, βdg, βmin, βmax)
    βc              = beta_at_equal_heights(a, S, V, βdg, βmin, βmax)
    @show run
    @show β0, βdg
    @show βl, βu
    @show βc
end