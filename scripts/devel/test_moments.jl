using LLRParsing
using Plots
using HDF5
using BenchmarkTools
using Profile
using LaTeXStrings
gr(size=(425,282),fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=10,legendfontsize=7,tickfontsize=7,labelfontsize=10,left_margin=0Plots.mm)

h5file = "data_assets/Sp4_Nt4_sorted.hdf5"
fid = h5open(h5file) 
r = keys(fid)[1]

a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid,r)
mina, maxa = extrema(a)
βs = collect(range(start=7.3402,stop=7.3405,length=100))

println("Timings for T = Float64")
@btime begin 
    LLRParsing.energy_moment.(Ref(S),Ref(a),βs,1,Float64)/(6V)^1;
    LLRParsing.energy_moment.(Ref(S),Ref(a),βs,2,Float64)/(6V)^2;
end
println("Timings for T = BigFloat (precision = $(precision(BigFloat(0))))")
@btime begin 
    LLRParsing.energy_moment.(Ref(S),Ref(a),βs,1,BigFloat)/(6V)^1;
    LLRParsing.energy_moment.(Ref(S),Ref(a),βs,2,BigFloat)/(6V)^2;
end

plt = plot()
for T in [Float64, BigFloat]
    E1 = LLRParsing.energy_moment.(Ref(S),Ref(a),βs,1,T)/(6V)^1;
    E2 = LLRParsing.energy_moment.(Ref(S),Ref(a),βs,2,T)/(6V)^2;
    CV = @. 6V*(E2 - E1^2)
    plot!(plt, βs, CV, label="$T", xlabel=L"\beta", ylabel=L"specific heat $C_V(\beta)$")
end
plt