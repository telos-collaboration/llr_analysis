using LLRParsing
using ArgParse
function lines_tex_body(file)
    lines = readlines(file)
    pattern = contains(raw"\hline \hline")
    l1 = findfirst(pattern, lines) + 1
    l2 = findlast(pattern, lines) - 1
    return lines[l1:l2]
end
function lines_tex_header(file)
    lines = readlines(file)
    pattern_header = contains(raw"\begin{tabular}")
    pattern_body = contains(raw"\hline \hline")
    l1 = findfirst(pattern_header, lines)
    l2 = findfirst(pattern_body, lines)
    return lines[l1:l2]
end
function lines_tex_footer(file)
    lines = readlines(file)
    pattern = contains(raw"\hline \hline")
    l1 = findlast(pattern, lines)
    return lines[l1:end]
end
function merge_tex_files(file1, file2, file_out)
    @assert lines_tex_header(file1) == lines_tex_header(file2)
    @assert lines_tex_footer(file1) == lines_tex_footer(file2)
    header = lines_tex_header(file1)
    footer = lines_tex_footer(file1)
    body1 = lines_tex_body(file1)
    body2 = lines_tex_body(file2)
    io = open(file_out, "w")
    print_provenance_tex(io)
    write.(Ref(io), header .* "\n")
    write.(Ref(io), body1 .* "\n")
    write.(Ref(io), "\\hline \\hline \n")
    write.(Ref(io), body2 .* "\n")
    write.(Ref(io), footer .* "\n")
    return close(io)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--file1"
        help = "First TeX file to be used"
        required = true
        "--file2"
        help = "Second TeX file to be used"
        required = true
        "--outfile"
        help = "Where to save the combined file"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file1 = args["file1"]
    file2 = args["file2"]
    outfile = args["outfile"]
    return merge_tex_files(file1, file2, outfile)
end
main()
