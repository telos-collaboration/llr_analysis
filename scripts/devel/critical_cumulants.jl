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

function cumulant_plots(h5file, Nt)
    fid = h5open(h5file)
    pltCV = plot(legend = :outerright, xlabel = L"\beta", ylabel = L"C_V(\beta)", title = L"specific heat $N_t = %$Nt$")
    pltBC = plot(legend = :outerright, xlabel = L"\beta", ylabel = L"B_L(\beta)", title = L"Binder cumulant $N_t = %$Nt$")
    pltCV_deriv = plot(legend = :outerright, xlabel = L"\beta", ylabel = L"C_V(\beta)", title = L"Derivative of specific heat $N_t = %$Nt$")
    pltBC_deriv = plot(legend = :outerright, xlabel = L"\beta", ylabel = L"B_L(\beta)", title = L"Derivative of Binder cumulant $N_t = %$Nt$")
    runs = largets_replica_runs(fid, keys(fid))

    for r in [last(runs)]

        a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid, r)
        β = range(start = minimum(a), stop = maximum(a), length = 200)
        repeats = size(a, 2)

        CV0 = zeros((length(β), repeats))
        BC0 = zeros((length(β), repeats))
        CV0_deriv = zeros((length(β), repeats))
        BC0_deriv = zeros((length(β), repeats))

        for i in 1:repeats
            EN(β, N) = LLRParsing.energy_moment(S, a[:, i], β, N, Float64, BigFloat) / (6.0V)^N
            E1(β) = EN(β, 1)
            E2(β) = EN(β, 2)
            E4(β) = EN(β, 4)
            fCV(β) = 6V * (E2(β) - E1(β)^2)
            fBC(β) = 1 - E4(β) / E2(β) / E2(β) / 3
            fCV_deriv(β, h = 0.001) = (- fCV(β + 2h) + 8fCV(β + h) - 8fCV(β - h) + fCV(β - 2h)) / (12h)
            fBC_deriv(β, h = 0.001) = (- fBC(β + 2h) + 8fBC(β + h) - 8fBC(β - h) + fBC(β - 2h)) / (12h)
            CV0[:, i] .= fCV.(β)
            BC0[:, i] .= fBC.(β)
            CV0_deriv[:, i] .= fCV_deriv.(β)
            BC0_deriv[:, i] .= fBC_deriv.(β)
        end

        CV = dropdims(mean(CV0, dims = 2), dims = 2)
        ΔCV = dropdims(std(CV0, dims = 2), dims = 2) ./ sqrt.(repeats)
        BC = dropdims(mean(BC0, dims = 2), dims = 2)
        ΔBC = dropdims(std(BC0, dims = 2), dims = 2) ./ sqrt.(repeats)

        CV_deriv = dropdims(mean(CV0_deriv, dims = 2), dims = 2)
        ΔCV_deriv = dropdims(std(CV0_deriv, dims = 2), dims = 2) ./ sqrt.(repeats)
        BC_deriv = dropdims(mean(BC0_deriv, dims = 2), dims = 2)
        ΔBC_deriv = dropdims(std(BC0_deriv, dims = 2), dims = 2) ./ sqrt.(repeats)

        plot!(pltCV, β, CV, ribbon = ΔCV, label = LLRParsing.fancy_title(r), lw = 2)
        plot!(pltBC, β, BC, ribbon = ΔBC, label = LLRParsing.fancy_title(r), lw = 2)
        plot!(pltCV_deriv, β, CV_deriv, ribbon = ΔCV_deriv, label = LLRParsing.fancy_title(r), lw = 2)
        plot!(pltBC_deriv, β, BC_deriv, ribbon = ΔBC_deriv, label = LLRParsing.fancy_title(r), lw = 2)
    end
    return pltCV, pltBC, pltCV_deriv, pltBC_deriv
end

Nt = 5
h5 = "data_assets/Sp4_Nt5_sorted.hdf5"

using BenchmarkTools

pltCV, pltBC, pltCV_deriv, pltBC_deriv = @time cumulant_plots(h5, Nt)
display(pltBC)
display(pltBC_deriv)
