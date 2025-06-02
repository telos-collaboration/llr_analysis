using HDF5
using LLRParsing
using Plots
using Statistics
using LaTeXStrings
using Peaks
using Roots
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

file  = "data_assets/test_Nt5_sorted.hdf5"
fid   = h5open(file)
betas = [7.48967, 7.48970, 7.48982, 7.48969, 7.48975]
runs  = keys(fid)
ind   = 5
beta  = betas[ind]
run   = runs[ind] 

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
    diff = peak_height_difference(fid,run,βold;w)

    for _ in 1:Nint
        βnew = diff>0 ? βold - βint/50 : βold + βint/50 
        newdiff = peak_height_difference(fid,run,βold;w)
        @show diff, newdiff
        @show βold, βnew
        # use extrame function to sort bracket of equal height
        if sign(diff) != sign(newdiff)
            return extrema((βold,βnew))
        end
        # move to the next interval
        βold, diff = βnew, newdiff
    end
end

using Roots
βl, βu = bracket_critical_beta(fid,run)

βl, βu = 7.48983354367, 7.489727429839999
peak_height_difference(fid,run,βl;w=50)
peak_height_difference(fid,run,βu;w=50)

find_zero(β -> peak_height_difference(fid,run,β), (βl, βu), Bisection())