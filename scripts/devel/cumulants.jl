using LLRParsing
using LaTeXStrings
using HDF5
using Plots
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

h5file = "data_assets/Sp4_Nt4_sorted.hdf5"
fid = h5open(h5file)
r = keys(fid)[end]

a, S, Nt, Nl, V = LLRParsing._set_up_histogram(fid, r)

repeats = size(a, 2)

aR = a[:, last(repeats)]
plaquette_moment(β, N) = LLRParsing.energy_moment(S, aR, β, N, BigFloat)


β = range(start = minimum(aR), stop = maximum(aR), length = 1000)
E1 = plaquette_moment.(β, 1) ./ ((6V)^1)
E2 = plaquette_moment.(β, 2) ./ ((6V)^2)
E4 = plaquette_moment.(β, 4) ./ ((6V)^4)

CV = @. 6V*(E2 - E1^2)
BC = @. 1 - E4/(3*E2^2)

pltCV = plot(xlabel = L"\beta", ylabel = "specirfic heat")
scatter!(pltCV, β, CV, label = L"C_V(\beta)")
pltBC = plot(xlabel = L"\beta", ylabel = "cumulants")
scatter!(pltBC, β[2:(end-1)], BC[2:(end-1)], label = L"C_B(\beta)")
