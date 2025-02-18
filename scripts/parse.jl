using BenchmarkTools
using HiRepParsing
using Parsers
function _parse_data!(array,string;n=6) 
    opts = Parsers.Options(delim=' ', ignorerepeated=true)
    io = IOBuffer(string)
    for i in 1:n
        array[i] = Parsers.parse(Float64, io, opts)
    end
end
function polyakov_loop(file)
    polyakov_loop = ComplexF64[]
    tmp = zeros(ComplexF32,2)
    pattern = "[FUND_POLYAKOV][0]Polyakov direction 0 = "
    i = length(pattern)
    for line in eachline(file)
        if startswith(line,pattern)
            _parse_data!(tmp,line[i:end];n=2)
            append!(polyakov_loop,tmp[1] + im*tmp[2])
        end
    end
    return polyakov_loop
end

file = "/home/fabian/Documents/Physics/Data/DataLLR/ImportanceSampling/Importance_Sampling_noCSV_noRW/IS_4x20/7.32/output_file"
T, L = latticesize(file)[1:2]
plaq = plaquettes(file)
beta = inverse_coupling(file)

poly = polyakov_loop(file)
poly = @profview polyakov_loop(file)
@btime polyakov_loop(file)
@show hash(poly)