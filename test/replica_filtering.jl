fid = h5open("./test_data/LLR_5x64_95_Rep13.hdf5")
fid_repeat = fid["5x64_1repeats_95replicas/13/"]

S0, steps, inds = remove_non_matching_trajectories_in_replicas(fid_repeat)
S0_sorted = sort(S0,dims=1)

Smin,Smax = extrema(filter(isfinite, S0))
nreplicas = first(size(S0))
S0_theory = collect(range(Smin,Smax,length=nreplicas))

@testset "Removal of corrupted data for non-matching trajectory lengths across replicas" begin 
    @test all(allequal ,eachslice(S0_sorted,dims=1))
    @test all(issorted ,eachslice(S0_sorted,dims=2))
    @test all(isapprox(S0_theory) ,eachslice(S0_sorted,dims=2))
end