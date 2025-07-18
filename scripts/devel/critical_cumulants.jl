using LLRParsing
using LaTeXStrings
using HDF5
using Plots
using ArgParse
using Statistics
using Quadmath

gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
)

function largets_replica_runs(h5id, runs)
    # Only include one run per volume with the largest number of N_replicas
    data = [[read(h5id[r], "Nt"), read(h5id[r], "Nl"), read(h5id[r], "N_replicas")] for r in runs]
    maxr = similar(runs)
    for i in eachindex(data)
        matches = findall(x -> x[1:2] == data[i][1:2], data)
        j = findmax(x -> data[x][3], matches)[2]
        maxr[i] = runs[matches[j]]
    end
    return unique(maxr)
end
function cumulants(h5dset, run)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(h5dset, run)
    β = range(start = minimum(a), stop = maximum(a), length = 100)
    return cumulants(h5dset, run, β)
end
function cumulants(h5dset, run, β)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(h5dset, run)
    repeats = size(a, 2)

    CV = zeros((length(β), repeats))
    BC = zeros((length(β), repeats))

    Threads.@threads for i in 1:repeats
        EN(β, N, V) = LLRParsing.energy_moment(S, a[:, i], β, N, Float128, BigFloat) / (6.0V)^N
        fCV(β) = 6V * (EN(β, 2, V) - EN(β, 1, V)^2)
        fBC(β) = 1 - EN(β, 4, 1) / EN(β, 2, 1) / EN(β, 2, 1) / 3
        CV[:, i] .= fCV.(β)
        BC[:, i] .= fBC.(β)
    end
    return β, CV, BC
end
function beta_extremal(β, obs; f = findmax)
    repeats = size(obs, 2)
    βmax0 = zeros(repeats)
    for i in 1:repeats
        ind = f(obs[:, i])[2]
        βmax0[i] = β[ind]
    end
    βmax = mean(βmax0)
    Δβmax = std(βmax0) / sqrt(repeats)
    return βmax, Δβmax
end

Nt = 4
setprecision(BigFloat, 106)
h5 = "data_assets/Sp4_Nt4_sorted.hdf5"

h5dset = h5open(h5)
runs = keys(h5dset)

function critical_beta_cumulants(h5dset, r; N = 200, eps = 1.0e-6, min_iter = 5, max_iter = 20, w = 20)
    a = first(LLRParsing._set_up_histogram(h5dset, r))
    min_a, max_a = minimum(a), maximum(a)
    β = range(start = min_a, stop = max_a, length = N)

    βc_old, Δβc_old = +Inf, +Inf
    for i in 1:max_iter
        β, CV0, BC0 = cumulants(h5dset, r, β)
        βc, Δβc = beta_extremal(β, CV0; f = findmax)
        βmin, βmax = extrema(β)
        βmin, βmax = min(βmin, βc - w * Δβc), max(βmin, βc + w * Δβc)
        β = range(start = (βmin + βc) / 2, stop = (βmax + βc) / 2, length = N)
        diff = abs(βc_old - βc)
        βc_old = βc
        if diff < eps && i > min_iter
            return βc, Δβc
        end
    end
    return βc, Δβc
end
