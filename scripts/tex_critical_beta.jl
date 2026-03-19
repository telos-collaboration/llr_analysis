using LLRParsing
using LaTeXStrings
using DelimitedFiles
using ArgParse
using Plots
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topright,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 1Plots.mm,
    palette = :Set1_5,
)
LINESTYLES = [:solid, :dot, :dash]
MARKERS = [:circle, :diamond, :dtriangle, :heptagon, :hexagon, :ltriangle, :octagon, :pentagon, :rect, :rtriangle, :star4, :star5, :star6, :star7, :star8, :utriangle ]

function read_critical_betas(file; offset = 0)
    data, header = readdlm(file, ',', header = true, comments = true)
    reps = Int.(data[:, 3])
    runs = data[:, 4]
    Nt = Int.(data[:, 5])
    Ns = Int.(data[:, 6])
    βc = Float64.(data[:, 9 + offset])
    Δβc = Float64.(data[:, 10 + offset])
    return runs, Nt, Ns, βc, Δβc, reps
end
function main_table(file1, file2, file_tex)
    runs11, Nt11, Ns11, βc11, Δβc11, reps11 = read_critical_betas(file1)
    runsCV, NtCV, NsCV, βcCV, ΔβcCV, repsCV = read_critical_betas(file2; offset = -2)
    runsBC, NtBC, NsBC, βcBC, ΔβcBC, repsBC = read_critical_betas(file2)

    io = open(file_tex, "w")
    print_provenance_tex(io)
    header = L"""\begin{tabular}{|c|c|c|c|c|c|} \hline
    $N_t$ & $N_s$ & $N_{\rm rep}$ & $\beta_{CV}(P)$ & $\beta_{CV}(C_V)$ & $\beta_{CV}(B_V)$ \\ \hline \hline"""
    footer = """\\hline \\hline
    \\end{tabular}"""

    @assert Nt11 == NtCV == NtBC
    println(io, header)
    for i in eachindex(Nt11, NtCV, NtBC)
        rx = r"[0-9]x[0-9]+_([0-9]+)replicas"
        m = match(rx, runs11[i])
        replicas = m.captures[1]
        str_11 = errorstring(βc11[i], Δβc11[i])
        str_CV = errorstring(βcCV[i], ΔβcCV[i])
        str_BC = errorstring(βcBC[i], ΔβcBC[i])
        Nt, Ns = Nt11[i], Ns11[i]
        println(io, "$Nt & $Ns & $replicas & $str_11 & $str_CV & $str_BC \\\\")
    end
    println(io, footer)
    return close(io)
end
function plot_critical_beta!(plt, Ns, βc, Δβc; xoffset, kws...)
    tks = (inv.(Ns), (L"1/%$Li" for Li in Ns))
    scatter!(plt, inv.(Ns) .+ xoffset, βc, xticks = tks, yerr = Δβc; kws...)
    return nothing
end
function largest_replica_number(Ns,reps)
    # filter out only runs with the smallest intervall
    inds = zeros(Int,length(unique(Ns)))
    
    for (i,L) in enumerate(unique(Ns))
        matches = findall(isequal(L),Ns)
        inds[i] = matches[findmax(reps[matches])[2]]
    end
    return inds
end
function main_plot(file1, file2, outfile)
    runs11, Nt11, Ns11, βc11, Δβc11, reps11 = read_critical_betas(file1)
    runsCV, NtCV, NsCV, βcCV, ΔβcCV, repsCV = read_critical_betas(file2; offset = -2)
    runsBC, NtBC, NsBC, βcBC, ΔβcBC, repsBC = read_critical_betas(file2)

    @assert Nt11 == NtCV == NtBC
    Nt = first(Nt11)
    plt = plot(title = L"N_t = %$Nt", ylabel = L"$\beta_{CV}$", xlabel = L"$1/N_s$", legend = :topleft)

    # only plot runs with fines intervals
    ind11 = largest_replica_number(Ns11,reps11)
    indCV = largest_replica_number(NsCV,repsCV)
    indBC = largest_replica_number(NsBC,repsBC)

    Δ = maximum(inv.(Ns11)) - minimum(inv.(Ns11))
    plot_critical_beta!(plt, NsBC[indBC], βcBC[indBC], ΔβcBC[indBC]; xoffset = 0 ,ma = 0.7, msa = 0.7, marker = :rect, label = L"$\beta_{CV }(B_C)$")
    plot_critical_beta!(plt, Ns11[ind11], βc11[ind11], Δβc11[ind11]; xoffset = -Δ/80 ,ma = 0.7, msa = 0.7, marker = :circ, label = L"$\beta_{CV }(P)$")
    plot_critical_beta!(plt, NsCV[indCV], βcCV[indCV], ΔβcCV[indCV]; xoffset = +Δ/80 ,ma = 0.7, msa = 0.7, marker = :diamond, label = L"$\beta_{CV }(C_V)$")
    savefig(outfile)

    return nothing
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--input_histogram"
        help = "CSV file with critical beta from histogram"
        required = true
        "--input_cumulants"
        help = "CSV file with critical beta from cumulants"
        required = true
        "--tex_file"
        help = "Where to save the table"
        required = true
        "--plot_file"
        help = "Where to save the table"
        required = true
    end
    return parse_args(s)
end
args = parse_commandline()
file1 = args["input_histogram"]
file2 = args["input_cumulants"]
file_tex = args["tex_file"]
file_plot = args["plot_file"]
main_table(file1, file2, file_tex)
main_plot(file1, file2, file_plot)
