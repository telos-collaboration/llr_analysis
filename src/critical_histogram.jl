"""
    peak_height_difference(args;w=5,nbins=length(S),A1=1,A2=1)

    Calculate the height difference between the two peaks in the histogram.

    `w` is the minimal number of intervals that the peaks need to be separated
    to count as proper peaks. This variable probably needs to be tuned depending
    on the noise in the dataset.

    The histogram is calculated at `nbins` values of the plaquette. 

    `A1` and `A2` set the relative heights of the peaks. They can
    be used to investigate histograms of unequal heights.
"""
function peak_height_difference(fid,run,β; kws...)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,run)
    return peak_height_difference(a, S, β, V; kws...)
end
function peak_height_difference(a, S, β, V; w=5, nbins=length(S),A1=1,A2=1,retries=0)
    ups, P, ΔP, covP, V, dS = probability_density(a, S, β, V; nbins)
    pks = findmaxima(P,w)
    n_peaks = length(pks.indices)
    if n_peaks == 2
        heightdiff = A2*pks.heights[2] - A1*pks.heights[1]
        return heightdiff
    elseif n_peaks > 2
        @warn "Found $n_peaks peak(s): Taking first and last one to be the relevant ones"
        heightdiff = A2*pks.heights[end] - A1*pks.heights[1]
        return heightdiff
    else
        @error "Found $n_peaks peak(s): Aborting."
    end
end
function _count_peaks(a, S, β, V; w=5, nbins=length(S), kws...)
    ups, P, ΔP, covP, V, dS = probability_density(a, S, β, V; nbins)
    pks = findmaxima(P,w)
    n_peaks = length(pks.indices)
    return n_peaks
end
"""
    find_initial_double_peak(a,S,β0,βmin,βmax,V; kwargs...)

    Find a value of β for which the histogram as tweo distinct peaks. 
    The starting value is β0 and we restrict ourselves to β ∈ [βmin,βmax]
"""
function find_initial_double_peak(a,S,β0,βmin,βmax,V;Nint=50, kws...)
    # start looking at β0 in the interval [βmin,βmax]
    np = _count_peaks(a, S, β0, V; kws...)
    if np == 2
        return β0
    end
    # If there is no double peak at β0, then partition the interval
    # into Nint intervals and check them in order of proximity to βmin.
    βs          = collect(range(βmin,βmax,length=Nint))
    ordered_ind = sort(eachindex(βs), by = ind -> abs(βs[ind] - β0))
    for i in ordered_ind
        np = _count_peaks(a, S, βs[i], V; kws...)
        if np == 2
            return βs[i]
        end
    end
end
"""
    bracket_critical_beta(a, S, V, β0, βmin, βmax; w=5, Nint=50)

    Find two value of beta for which the height difference has a different sign.
    The two values can then be used to determine the critical beta using a Bisection
    search to find the value for which the height difference is vanishing.

    `w` is the minimal required separation of the peaks.

    Here I start with an initial β0 for which we have to peaks. Then, depending on 
    the sign of the height difference, β is lowered or raised in fixed steps. The maximal 
    number of steps is given by `Nint`.    
"""
function bracket_critical_beta(a, S, V, β0, βmin, βmax;Nint=100, kws...)
    βold = β0
    βint = abs(βmax - βmin)/2
    for _ in 1:Nint
        diff    = peak_height_difference(a, S, βold, V; kws... )
        βnew    = diff>0 ? βold - βint/Nint : βold + βint/Nint 
        newdiff = peak_height_difference(a, S, βnew, V; kws... )
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
function beta_at_equal_heights(fid,run,β0, βmin, βmax;kws...)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,run)
    return beta_at_equal_heights(a, S, V, β0, βmin, βmax;kws...)
end
function beta_at_equal_heights(a, S, V;kws...)
    βmin, βmax      = extrema(a)
    β0              = (βmin + βmax)/2
    return beta_at_equal_heights(a, S, V, β0, βmin, βmax;kws...)
end
function beta_at_equal_heights(a, S, V, β0, βmin, βmax; kws...)
    βdg    = find_initial_double_peak(a,S,β0,βmin,βmax,V; kws...)
    βl, βu = bracket_critical_beta(a, S, V, βdg, βmin, βmax;kws...)
    f(β)   = peak_height_difference(a, S, β, V;kws...)
    βc     = find_zero(f, (βl, βu), Bisection())
    return βc
end