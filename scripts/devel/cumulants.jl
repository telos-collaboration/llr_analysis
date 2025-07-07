# In python, we are using the following precision
# Mpmath settings:
#   mp.prec = 53                [default: 53]
#   mp.dps = 15                 [default: 15]
#   mp.trap_complex = False     [default: False]
setprecision(BigFloat, 53)
using LLRParsing

function log_partition_function(a, S, beta)
    # David uses a different sign for a
    dS     = S[2] - S[1]    
    pi_exp = 0.0
    Z      = BigFloat(0)
    for (ai, Si) in zip(a,S)
        A           = beta - ai
        exp_factor  = exp(BigFloat(pi_exp + Si*beta - ai*dS/2))
        sinh_factor = iszero(A) ? dS/2 : sinh(A*dS/2)/A
        Z      += exp_factor*sinh_factor
        pi_exp -= ai*dS
    end
    return log(2Z)
end

"""
    Calculate n-th cumulant of the energy from the d.o.s.
    This code has been ported over from David Mason's implementation
    in python (10.5281/zenodo.13807993). It follows Eqs. (3.1.18) and
    (3.1.19) from David Masons's PhD thesis. 
"""
function energy_cumulant(S,a,β,N::Int)
    pi_exp = - log_partition_function(a, S, β)
    full_exp = BigFloat(0)
    En = BigFloat(0)
    δS = S[2] - S[1]
    for (Si,ai) in zip(S,a)
        A = - BigFloat(ai) + BigFloat(β)
        
        full_exp = 2*factorial(N)*exp(pi_exp + β*(Si-δS/2) + A*δS/2 ) 
        T = A*δS/2
        for m in 0:N
            sinh_term = BigFloat(0)
            cosh_term = BigFloat(0)
            for j in 0:div(m,2,RoundDown)
                sinh_term += (δS/2)^(2j)*Si^(m-2j)/factorial(2j)/factorial(m-2j)
            end
            for j in 1:div(m,2,RoundUp)
                cosh_term += (δS/2)^(2j-1)*Si^(m-2j+1)/factorial(2j-1)/factorial(m-2j+1)
            end
            En += full_exp*(sinh(T)*sinh_term + cosh(T)*cosh_term)*A^(m-N-1)*(-1)^(N-m)
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

h5file = "data_assets/Sp4_Nt4_sorted.hdf5"
fid = h5open(h5file) 
r = keys(fid)[1]

a, S, Nt, Nl, V = _set_up_histogram(fid,r)
mina, maxa = extrema(a)
βs = range(start=mina,stop=maxa,length=100)

E0 = [ energy_cumulant(S,a,β,0)/(6V)^0 for β in βs ]
E1 = [ energy_cumulant(S,a,β,1)/(6V)^1 for β in βs ]
E2 = [ energy_cumulant(S,a,β,2)/(6V)^2 for β in βs ]
E4 = [ energy_cumulant(S,a,β,4)/(6V)^4 for β in βs ]

CV = @. 6V*(E2 - E1^2)
plot(βs, CV)
