using LLRParsing
using LaTeXStrings
using HDF5
using Plots
using ArgParse

using Statistics
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

function main(h5file, β)
    fid = h5open(h5file)
    pltCV = plot(legend = :outerright, xlabel = L"\beta", ylabel = L"C_V(\beta)", title = "specific heat")
    pltBC = plot(legend = :outerright, xlabel = L"\beta", ylabel = L"B_L(\beta)", title = "Binder cumulant")
    runs = largets_replica_runs(fid, keys(fid))

    for r in runs

        a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid, r)
        repeats = size(a, 2)

        CV0 = zeros((length(β), repeats))
        BC0 = zeros((length(β), repeats))

        for i in 1:repeats
            E1 = LLRParsing.energy_moment.(Ref(S), Ref(a[:, i]), β, 1) ./ ((6V)^1)
            E2 = LLRParsing.energy_moment.(Ref(S), Ref(a[:, i]), β, 2) ./ ((6V)^2)
            E4 = LLRParsing.energy_moment.(Ref(S), Ref(a[:, i]), β, 4) ./ ((6V)^4)
            @. CV0[:, i] = 6V * (E2 - E1^2)
            @. BC0[:, i] = 1 - E4 / (3 * E2^2)
        end


        CV = dropdims(mean(CV0, dims = 2), dims = 2)
        ΔCV = dropdims(std(CV0, dims = 2), dims = 2) ./ sqrt.(repeats)
        BC = dropdims(mean(BC0, dims = 2), dims = 2)
        ΔBC = dropdims(std(BC0, dims = 2), dims = 2) ./ sqrt.(repeats)

        plot!(pltCV, β, CV, ribbon = ΔCV, label = LLRParsing.fancy_title(r), lw = 2)
        plot!(pltBC, β[2:(end - 1)], BC[2:(end - 1)], ribbon = ΔBC[2:(end - 1)], label = LLRParsing.fancy_title(r), lw = 2)
        display(pltCV)
        display(pltBC)
    end
    return pltCV, pltBC
end

h5file = "data_assets/Sp4_Nt4_sorted.hdf5"
β = range(start = 7.338, stop = 7.342, length = 100)
pltCV, pltBC = main(h5file, β)
#savefig(pltCV,"specific_heat_Nt4.pdf")

#h5file = "data_assets/Sp4_Nt5_sorted.hdf5"
#β      = range(start=7.489,stop=7.4905,length=100)
#pltCV, pltBC  = main(h5file,β)
##savefig(pltCV,"specific_heat_Nt5.pdf")
