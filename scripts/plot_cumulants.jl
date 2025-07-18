using LLRParsing
using LaTeXStrings
using HDF5
using Plots
using ArgParse
using Statistics
using Peaks
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
function cumulant_plots(h5file, Nt)
    fid = h5open(h5file)
    pltCV = plot(legend = :outerright, xlabel = L"\beta", ylabel = L"C_V(\beta)", title = L"specific heat $N_t = %$Nt$")
    pltBC = plot(legend = :outerright, xlabel = L"\beta", ylabel = L"B_L(\beta)", title = L"Binder cumulant $N_t = %$Nt$")
    runs = largets_replica_runs(fid, keys(fid))

    for r in runs
        β, CV0, BC0 = cumulants(fid, r)
        repeats = size(CV0, 2)

        CV = dropdims(mean(CV0, dims = 2), dims = 2)
        ΔCV = dropdims(std(CV0, dims = 2), dims = 2) ./ sqrt.(repeats)
        BC = dropdims(mean(BC0, dims = 2), dims = 2)
        ΔBC = dropdims(std(BC0, dims = 2), dims = 2) ./ sqrt.(repeats)

        plot!(pltCV, β, CV, ribbon = ΔCV, label = LLRParsing.fancy_title(r), lw = 2)
        plot!(pltBC, β, BC, ribbon = ΔBC, label = LLRParsing.fancy_title(r), lw = 2)
    end
    hline!(pltBC, [2 / 3], yformatter = :plain, linestyle = :dot, color = :black, label = L"B_L(\beta) = 2/3")
    return pltCV, pltBC
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--plot_file_specific_heat"
        help = "Where to save the plot"
        required = true
        "--plot_file_binder_cumulant"
        help = "Where to save the plot"
        required = true
        "--Nt"
        help = "Nt of the runs to be plotted of the plot"
        required = true
        arg_type = Int
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    pltCV, pltBC = cumulant_plots(args["h5file"], args["Nt"])
    savefig(pltCV, args["plot_file_specific_heat"])
    savefig(pltBC, args["plot_file_binder_cumulant"])
    return nothing
end
main()
