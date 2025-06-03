using HDF5
using LLRParsing
using Plots
using Statistics
using LaTeXStrings
using Peaks
using Roots
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function peak_height_difference(fid,run,β;w=50)
    P   = probability_density(fid, run, β)[2]
    pks = findmaxima(P,w)
    @assert length(pks.indices) == 2
    heightdiff = pks.heights[2] - pks.heights[1]
    return heightdiff
end
function bracket_critical_beta(fid,run;w=50,Nint=50)
    a = LLRParsing.a_vs_central_action_repeats(fid,run;ind=nothing)[1]
    βl, βu = extrema(a)
    βold = (βl + βu)/2
    βint = abs(βl - βu)/2
    for _ in 1:Nint
        diff    = peak_height_difference(fid,run,βold;w)
        βnew    = diff>0 ? βold - βint/Nint : βold + βint/Nint 
        newdiff = peak_height_difference(fid,run,βnew;w)

        if sign(diff) != sign(newdiff)
            # use extrame function to sort bracket of equal height
            return extrema((βold,βnew))
        end
        # move to the next interval
        βold, diff = βnew, newdiff
    end
end
function beta_at_equal_heights(fid,run)
    βl, βu = bracket_critical_beta(fid,run)
    βc = find_zero(β -> peak_height_difference(fid,run,β), (βl, βu), Bisection())
    return βc
end

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
runs  = keys(fid)
βc    = @profview beta_at_equal_heights(fid,last(runs))