using LibGit2
using Dates
using HDF5
using InteractiveUtils
function provenance()
    provenance = Dict{String, String}()

    provenance["comment"] = """
    This file and all the files in this directory were generated automatically.
    Do not modify them; re-run the analysis workflow!"""

    # filename and command-line parameters (includes all input)
    provenance["script"] = @__FILE__
    provenance["cli_args"] = prod(ARGS .* " ")

    # git has and status of the repo
    try
        githash = LibGit2.head(".")
        githash = LibGit2.isdirty(LibGit2.GitRepo(".")) ? githash * "-dirty" : githash
        provenance["git_hash"] = githash
    catch
        provenance["git_hash"] = "[No commit ID available]"
    end

    # get versioninfo for julia and machine
    buf = IOBuffer()
    InteractiveUtils.versioninfo(buf)
    provenance["versioninfo"] = String(take!(buf))[1:(end - 1)]

    # current time, hostname and user
    provenance["timestamp"] = string(now())
    provenance["machine"] = gethostname()
    provenance["user"] = ENV["USER"]

    return provenance
end
function print_provenance_tex(io = stdout)
    comment_str = "% "
    return print_provenance_plaintext(io, comment_str)
end
function print_provenance_csv(io = stdout)
    comment_str = "# "
    return print_provenance_plaintext(io, comment_str)
end
function print_provenance_plaintext(io = stdout, comment_str = "% ")
    p = provenance()
    str = string(p["comment"], "\n") *
        string("workflow script: ", p["script"], "\n") *
        string("workflow input and args: ", p["cli_args"], "\n") *
        string("workflow git hash: ", p["git_hash"], "\n") *
        string("julia versioninfo: ", p["versioninfo"], "\n") *
        string("workflow timestamp: ", p["timestamp"], "\n") *
        string("analysis machine: ", p["machine"], "\n") *
        string("analysis user: ", p["user"], "\n")
    str = replace(comment_str * str, "\n" => "\n" * comment_str)[1:(end - 2)]
    return print(io, str)
end
function write_provenance_hdf5(h5file)
    p = provenance()
    for k in keys(p)
        h5write(h5file, joinpath("provenance", k), p[k])
    end
    return
end
