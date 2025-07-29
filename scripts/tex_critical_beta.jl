using LLRParsing
using LaTeXStrings
using DelimitedFiles
using ArgParse

function read_critical_betas(file; offset = 0)
    data = readdlm(file, ',', skipstart = 1)
    runs = data[:, 1]
    Nt = data[:, 2]
    Nl = data[:, 3]
    βc = data[:, 6 + offset]
    Δβc = data[:, 7 + offset]
    return runs, Nt, Nl, βc, Δβc
end
function main(file1, file2, file_tex)
    runs11, Nt11, Nl11, βc11, Δβc11 = read_critical_betas(file1)
    runsCV, NtCV, NlCV, βcCV, ΔβcCV = read_critical_betas(file2; offset = -2)
    runsBC, NtBC, NlBC, βcBC, ΔβcBC = read_critical_betas(file2)

    io = open(file_tex, "w")
    header = L"""\begin{tabular}{|c|c|c|c|c|c|} \hline
    $N_t$ & $N_l$ & $N_{\rm rep}$ & $\beta^{1:1}_c$ & $\beta^{C_V}_c$ & $\beta^{B_V}_c$ \\ \hline \hline """
    footer = """\\hline \\hline
    \\end{tabular}"""

    @assert Nt11 == NtCV == NtBC
    @assert Nl11 == NlCV == NlBC
    println(io, header)
    for i in eachindex(Nt11, NtCV, NtBC)
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
    return close(io)
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
    end
    return parse_args(s)
end
args = parse_commandline()
file1 = args["input_histogram"]
file2 = args["input_cumulants"]
file_tex = args["tex_file"]
main(file1, file2, file_tex)
