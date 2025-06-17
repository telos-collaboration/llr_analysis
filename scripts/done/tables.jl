using Pkg; Pkg.activate(".")
using LLRParsing
using HDF5
using LaTeXStrings
using ArgParse

function write_run_table(file,outfile)
    fid = h5open(file)
    io  = open(outfile,"w")

    header = L"\begin{tabular}{|c|c|c|c|c|c|c|c|} \hline 
$N_t$ & $N_l$ & $u_{p}^{\rm min}$ & $u_{p}^{\rm max}$ & $N_{\rm rep}$ & $N_{\rm repeats}$ & $n_{\rm NR}$ & $n_{\rm NR}$ \\ \hline \hline "
    footer = "\\hline \\hline 
\\end{tabular}"

    runs = keys(fid)
    println(io,header)
    for run in runs
        a0, Δa0, S0, ind = a_vs_central_action(fid,run)
        Nt = read(fid[run],"Nt")
        Nl = read(fid[run],"Nl")
        repeats   = read(fid[run],"repeats")
        Nrepeats  = read(fid[run],"N_repeats")
        Nreplicas = read(fid[run],"N_replicas")
        up = S0/(6Nl^3 * Nt)
        isrm = read(fid[run],"$(first(repeats))/Rep_0/is_rm")
        nr, rm = findlast(x->!x,isrm), length(isrm)
        println(io,"$Nt & $Nl & $(minimum(up)) & $(maximum(up)) & $Nreplicas & $Nrepeats & $nr & $(rm-nr) \\\\")
    end
    println(io,footer)
    close(io)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
            help = "HDF5 file containing the sorted results"
            required = true
        "--outfile"
            help = "Where to save the table"
            required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    write_run_table(args["h5file"],args["outfile"])
end
main()