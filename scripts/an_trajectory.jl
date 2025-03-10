function sort_by_central_energy_to_hdf5(h5file_in,h5file_out)
    h5dset = h5open(h5file_in)
    for run in keys(h5dset)
        N_replicas = read(h5dset[run],"N_replicas")
        N_repeats  = read(h5dset[run],"N_repeats")
        # read all last elements for a and the central action
        for j in 1:N_repeats
            ntraj = length(h5dset[run]["$(j-1)/Rep_0/is_rm"])
            a = zeros(N_replicas,ntraj)
            S = zeros(N_replicas,ntraj)
            p = zeros(N_replicas,ntraj)
            for i in 1:N_replicas
                a[i,:] = h5dset[run]["$(j-1)/Rep_$(i-1)/a"][] 
                S[i,:] = h5dset[run]["$(j-1)/Rep_$(i-1)/S0"][]
                p[i,:] = h5dset[run]["$(j-1)/Rep_$(i-1)/plaq"][]
            end
            ## Sort by the central action to account for different swaps
            for j in 1:ntraj
                perm = sortperm(S[:,j])
                S[:,j] = S[perm,j]
                a[:,j] = a[perm,j]
                p[:,j] = p[perm,j]
            end
            ## make sure that the sorted central action alwas matches
            for i in 1:N_replicas
                @assert allequal(S[i,:])
                h5write(h5file_out,joinpath(run,"$(j-1)","Rep_$(i-1)","S0_sorted"),  S[i,:])
                h5write(h5file_out,joinpath(run,"$(j-1)","Rep_$(i-1)","a_sorted"),   a[i,:])
                h5write(h5file_out,joinpath(run,"$(j-1)","Rep_$(i-1)","plaq_sorted"),p[i,:])
                h5write(h5file_out,joinpath(run,"$(j-1)","Rep_$(i-1)","dS0"),h5read(h5file_in,joinpath(run,"$(j-1)","Rep_$(i-1)","dS0")))
            end
        end
        h5write(h5file_out,joinpath(run,"N_replicas"),N_replicas)
        h5write(h5file_out,joinpath(run,"N_repeats"),N_repeats)
        h5write(h5file_out,joinpath(run,"Nt"),h5read(h5file_in,joinpath(run,"Nt")))
        h5write(h5file_out,joinpath(run,"Nl"),h5read(h5file_in,joinpath(run,"Nl")))
    end
end

using HDF5

h5file_in = "output/SU3_llr_david.hdf5"
h5file_out= "output/SU3_llr_david_sorted.hdf5"
sort_by_central_energy_to_hdf5(h5file_in,h5file_out)

#h5file_in = "output/Sp4_llr_david.hdf5"
#h5file_out= "output/Sp4_llr_david_sorted.hdf5"
#sort_by_central_energy_to_hdf5(h5file_in,h5file_out)