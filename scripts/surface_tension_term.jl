using LLRParsing
using HDF5
using Peaks
using Statistics
using Plots
using LaTeXStrings
using ArgParse
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

function beta_Pmin_Pmax_jackknife(fid, run; w = 5)
    a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid, run)
    a_jk = LLRParsing.jackknife_resamples(a)

    beta = zeros(size(a_jk, 2))
    Pmin = zeros(size(a_jk, 2))
    Pmax = zeros(size(a_jk, 2))
    inter = zeros(size(a_jk, 2))

    for i in axes(a_jk, 2)
        ai = a_jk[:, i:i]
        beta[i] = LLRParsing.beta_at_equal_heights(ai, S, V)
        ups, P, ΔP = probability_density(ai, S, beta[i], V)[1:3]
        # find the two peaks of the probability distribution
        # make sure that we have found two peaks
        # find the minimum in between the peaks
        pks = findmaxima(P, w)
        @assert length(pks.indices) == 2
        i1, i2 = pks.indices
        Pmin[i] = minimum(P[i1:i2])
        Pmax[i] = (P[i1] + P[i2]) / 2
        # interface term according to Eq.(28) in Bennett:2024bhy
        # Note that the second term has a typo: the sign is flipped
        x = (Nt / Nl)^2
        inter[i] = - x * log(Pmin[i] / Pmax[i]) / 2 + x * log(Nl) / 4
    end
    return beta, Pmin, Pmax, inter, Nt, Nl
end
function apply_jackknife(obs::AbstractVector)
    N = length(obs)
    O = mean(obs)
    ΔO = sqrt(N - 1) * std(obs, corrected = false)
    return O, ΔO
end
function main(files, plt_name)
    plt = plot(; ylabel = L"\tilde{I}", xlabel = L"N_t^2/N_s^2")
    Nt = 0
    for file in files
        fid = h5open(file)
        runs = keys(fid)
        runs = filter(!startswith("provenance"), runs)
        x, I, ΔI = zeros(length(runs)), zeros(length(runs)), zeros(length(runs))
        for (i, r) in enumerate(runs)
            beta, Pmin, Pmax, inter, Nt, Nl = beta_Pmin_Pmax_jackknife(fid, r)
            I[i], ΔI[i] = apply_jackknife(inter)
            x[i] = inv(Nl / Nt)
        end
        plot!(plt, x .^ 2, I, yerr = ΔI, markershape = :circle, markeralpha = 0.7, label = L"N_t=%$Nt")
    end
    plot!(plt; ylims = (0, maximum(ylims(plt))))
    plot!(plt; xlims = (0, maximum(xlims(plt))))
    savefig(plt, plt_name)
    return plt
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--plotfile"
        help = "Where to save the plot"
        required = true
        "arg"
        help = "HDF5 files with sorted data for all Nt to be plotted"
        required = true
        nargs = '+'
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    plotfile = args["plotfile"]
    files = args["arg"]
    main(files, plotfile)
    return nothing
end
main()
