using HDF5
using LLRParsing
using Plots
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

file = "data_assets/test_Nt5_sorted.hdf5"
fid  = h5open(file)
runs = keys(fid)

title = L"contination $5\times48^3$"
a0    = hcat([read(fid[runs[1]],"0/Rep_$i/a_sorted") for i in 0:47]...)
plt   = plot(a0;ms=1,legend=:outerright,title,label="",xlabel="updates (NR + RM)",ylabel=L"a_n")
savefig("Nt5_Nl48_continued.pdf")