# In python, we are using the following precision
# Mpmath settings:
#   mp.prec = 53                [default: 53]
#   mp.dps = 15                 [default: 15]
#   mp.trap_complex = False     [default: False]
setprecision(BigFloat, 53)
using LLRParsing

function log_partition_function(a, S, beta, ::Type{T}=BigFloat) where T
    # David uses a different sign for a
    dS     = S[2] - S[1]    
    pi_exp = T(0)
    Z      = T(0)
    for (ai, Si) in zip(a,S)
        A           = beta - ai
        exp_factor  = exp(T(pi_exp + Si*beta - ai*dS/2))
        sinh_factor = iszero(A) ? dS/2 : sinh(A*dS/2)/A
        Z      += exp_factor*sinh_factor
        pi_exp -= ai*dS
    end
    return Float64(log(2Z))
end

"""
    llr_energy_momenta(S,a,β,N,[::Type{T}])

    Calculate N-th cumulant of the energy from the d.o.s. coefficients `a`
    at the central energies `S` at the inverse coupling β.

    This code has been ported over from David Mason's implementation
    in python (10.5281/zenodo.13807993). It follows Eqs. (3.1.18) and
    (3.1.19) from David Masons's PhD thesis. 

    By default, we use julia's BigFloat together with the double precision Float64 
    datatype for floating point calculations. This can be overwritten by specifying
    the datatype `T` to be used instead.

    Typically, it suffices to use Float64 for all calculations in here as long 
    as the log of the partition function is calculated in higher precision.
"""
function llr_energy_momenta(S,a,β,N::Int,::Type{T}=Float64) where T
    pi_exp = - T(log_partition_function(a, S, β, BigFloat))
    full_exp = T(0)
    En = T(0)
    δS = S[2] - S[1]
    for (Si,ai) in zip(S,a)
        A = T(- ai + β)
        full_exp = exp(pi_exp + β*(Si-δS/2) + A*δS/2 ) 
        for m in 0:N
            sinh_term = T(0)
            cosh_term = T(0)
            for j in 0:div(m,2,RoundDown)
                sinh_term += (δS/2)^(2j)*Si^(m-2j)/factorial(2j)/factorial(m-2j)
            end
            for j in 1:div(m,2,RoundUp)
                cosh_term += (δS/2)^(2j-1)*Si^(m-2j+1)/factorial(2j-1)/factorial(m-2j+1)
            end
            sh = sinh(A*δS/2)
            ch = cosh(A*δS/2)
            Ap = A^(m-N-1)
            En += 2*factorial(N)*full_exp*(sh*sinh_term + ch*cosh_term)*Ap*(-1)^(N-m)
        end   
        pi_exp -= ai*δS
    end
    return En
end
function _set_up_histogram(fid,run)
    a, S_all = LLRParsing.a_vs_central_action_repeats(fid,run;ind=nothing)[1:2]
    Nl   = read(fid[run],"Nl")
    Nt   = read(fid[run],"Nt")
    S    = unique(S_all) 
    V    = Nt*Nl^3
    return a, S, Nt, Nl, V
end

using Plots
using HDF5
using BenchmarkTools
using Profile
using LaTeXStrings
gr(size=(425,282),fontfamily="Computer Modern",legend=:topleft,frame=:box,titlefontsize=10,legendfontsize=7,tickfontsize=7,labelfontsize=10,left_margin=0Plots.mm)

h5file = "data_assets/Sp4_Nt4_sorted.hdf5"
fid = h5open(h5file) 
r = keys(fid)[1]

a, S, Nt, Nl, V = _set_up_histogram(fid,r)
mina, maxa = extrema(a)
βs = collect(range(start=7.3402,stop=7.3405,length=1000))

println("Timings for T = Float64")
@btime begin 
    llr_energy_momenta.(Ref(S),Ref(a),βs,1,Float64)/(6V)^1;
    llr_energy_momenta.(Ref(S),Ref(a),βs,2,Float64)/(6V)^2;
end
println("Timings for T = BigFloat (precision = $(precision(BigFloat(0))))")
@btime begin 
    llr_energy_momenta.(Ref(S),Ref(a),βs,1,BigFloat)/(6V)^1;
    llr_energy_momenta.(Ref(S),Ref(a),βs,2,BigFloat)/(6V)^2;
end

plt = plot()
for T in [Float64, BigFloat]
    E1 = llr_energy_momenta.(Ref(S),Ref(a),βs,1,T)/(6V)^1;
    E2 = llr_energy_momenta.(Ref(S),Ref(a),βs,2,T)/(6V)^2;
    CV = @. 6V*(E2 - E1^2)
    plot!(plt, βs, CV, label="$T", xlabel=L"\beta", ylabel=L"specific heat $C_V(\beta)$")
end
plt