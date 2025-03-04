using Pkg; Pkg.activate(".")
using LLRParsing

h5file = "LLR_data.hdf5"
ens    = "ImportanceSampling/4x20/7.32/"
binder, Δbinder, sh_plaq, Δsh_plaq, plaq, Δplaq, poly_sus, Δpoly_sus, poly, Δpoly = std_observables(h5file,ens)