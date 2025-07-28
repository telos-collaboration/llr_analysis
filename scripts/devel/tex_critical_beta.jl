using LLRParsing
using LaTeXStrings
using Plots
using DelimitedFiles
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

file1 = "data_assets/critical_beta_Nt5.csv"
file2 = "data_assets/critical_beta_2:1_Nt5.csv"
file3 = "data_assets/critical_beta_cumulants_Nt5.csv"

function read_critical_betas(file; offset = 0)
    data = readdlm(file, ',', skipstart = 1)
    runs = data[:, 1]
    Nt = data[:, 2]
    Nl = data[:, 3]
    βc = data[:, 6 + offset]
    Δβc = data[:, 7 + offset]
    return runs, Nt, Nl, βc, Δβc
end

runs11, Nt11, Nl11, βc11, Δβc11 = read_critical_betas(file1)
runs21, Nt21, Nl21, βc21, Δβc21 = read_critical_betas(file2)
runsCV, NtCV, NlCV, βcCV, ΔβcCV = read_critical_betas(file3; offset = -2)
runsBC, NtBC, NlBC, βcBC, ΔβcBC = read_critical_betas(file3)

plt = plot()
scatter!(plt, inv.(Nl11 ./ Nt11), βc11, yerr = Δβc11, label = "")
scatter!(plt, inv.(NlCV ./ NtCV), βcCV, yerr = ΔβcCV, label = "")
scatter!(plt, inv.(NlBC ./ NtBC), βcBC, yerr = ΔβcBC, label = "")

io = stdout

header = L"\begin{tabular}{|c|c|c|c|c|c|} \hline
$N_t$ & $N_l$ & $N_{\rm rep}$ & $\beta^{1:1}_c$ & $\beta^{C_V}_c$ & $\beta^{B_V}_c$ \\ \hline \hline "
footer = "\\hline \\hline
\\end{tabular}"


@assert Nt11 == Nt21 == NtCV == NtBC
@assert Nl11 == Nl21 == NlCV == NlBC
println(io, header)
for i in eachindex(Nt11, Nt21, NtCV, NtBC)
    rx = r"[0-9]x[0-9]+_([0-9]+)replicas"
    m = match(rx, runs11[i])
    replicas = m.captures[1]
    str_11 = errorstring(βc11[i], Δβc11[i])
    str_CV = errorstring(βcCV[i], ΔβcCV[i])
    str_BC = errorstring(βcBC[i], ΔβcBC[i])
    Nt, Nl = Nt11[i], Nl11[i]
    println(io, "$Nt & $Nl & $replicas & $str_11 & $str_CV & $str_BC \\\\")
end
println(io, footer)
