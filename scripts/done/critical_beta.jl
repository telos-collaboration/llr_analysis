using LLRParsing
using HDF5
using Plots
using ArgParse
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function all_critical_beta(file, outfile, β0, βmin, βmax; A1=1, A2=1)
    fid    = h5open(file)
    runs   = keys(fid)
    header = !isfile(outfile)
    io     = open(outfile,"a")
    header && println(io,"run,Nt,Nl,A1,A2,βc,Δβc,str")
    for run in runs
        L   = read(fid[run],"Nl")
        T   = read(fid[run],"Nt")
        try
            βc, Δβc = LLRParsing.βc_jackknife(fid,run;β0,βmin,βmax,A1,A2)
            str = errorstring(βc, Δβc; nsig=1)
            println(io,"$run,$T,$L,$A1,$A2,$βc,$Δβc,$str")
        catch 
            @warn "Cannot determin critical coupling for run $run"
            println(io,"$run,$T,$L,$A1,$A2,NaN,NaN,NaN(NaN)")
        end
    end
    close(io)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
            help = "HDF5 file containing the sorted results"
            required = true
        "--outfile"
            help = "Where to save results"
            required = true
        "--beta_min"
            help = "lower limit of β for finding the critical coupling"
            arg_type = Float64
            required = true
        "--beta_max"
            help = "upper limit of β for finding the critical coupling"
            arg_type = Float64
            required = true
        "--relative_peak_height"
            help = "Ratio of peak heights used for finding the critical coupling"
            arg_type = Int
            default = 1
        end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    βmin = args["beta_min"]
    βmax = args["beta_max"]
    file = args["h5file"]
    A1   = args["relative_peak_height"]
    β0   = (βmax+βmin)/2
    outfile = args["outfile"]
    isfile(outfile) && rm(outfile)
    all_critical_beta(file, outfile, β0, βmin, βmax; A1=A1, A2=1)
end
main()