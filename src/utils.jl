export combinePRS,
    combinePRSbyType,
    makeModelFrameCC,
    estimate_h2,
    estimate_neff,
    get_chi_markers,
    postfix_names!,
    snp_join,
    read_chi_markers

using DataFrames
using CSV


    

"add up PRS scores from chromosome files in dir and write to outfile"
function combinePRS(dir::AbstractString, outfile)

    columns = []
    fh = open(joinpath(dir, "prstypes.txt"))
    for line in eachline(fh)
        push!(columns, Symbol(chomp(line)))
    end
    close(fh)
    
    prs = Dict{String, Array{Float64}}()
    files = [joinpath(dir, f) for f in readdir(dir) if endswith(f, ".prs")]
    for f in files
        for line in eachline(f)
            cols = split(chomp(line), "\t")
            pn = cols[1]
            vals = parse.(Float64, cols[5:end])
            haskey(prs, pn) || (prs[pn] = zeros(Float64, length(vals)))
            prs[pn] .+= vals
        end
    end

    fh = open(outfile, "w")
    columns == [] || println(fh, join(vcat([:PN], columns), "\t"))
    for (k,vs) in sort(collect(prs))
        println(fh, k, "\t", join([string(round(v; sigdigits=4)) for v in vs], "\t"))
    end
    close(fh)
end


"for each type of prs weight in columns, create file containing weights by chromosome in columns"
function combinePRSbyType(dir::AbstractString, outdir)

    columns = []
    fh = open(joinpath(dir, "prstypes.txt"))
    for line in eachline(fh)
        push!(columns, Symbol(chomp(line)))
    end
    close(fh)

    chroms = ["chr" * string(i) for i = 1:22]
    
    prs = Dict{String, Dict{String, Array{Float64}}}()
    for chrom in chroms
        prs[chrom] = Dict{String, Array{Float64}}()  # prs[chrom][pn][type] = weight
        prs_chrom = prs[chrom]
        path = joinpath(dir, chrom * ".prs")
        for line in eachline(path)
            cols = split(chomp(line), "\t")
            pn = cols[1]
            vals = parse.(Float64, cols[5:end])
            @assert length(vals) == length(columns) "Wrong number of columns for $pn in $path"
            prs_chrom[pn] = vals
        end
    end

    idx = 0
    for column in columns
        idx += 1
        path = joinpath(outdir, string(column) * ".prs")
        fh = open(path, "w")
        println(fh, "PN\t", join([string(column) * "_" *string(chrom) for chrom in chroms], "\t"))
        pns = sort(collect(keys(prs[chroms[1]])))
        for pn in pns
            println(fh, pn, "\t", join([string(round(prs[chrom][pn][idx]; sigdigits=4)) for chrom in chroms], "\t"))
        end
        close(fh)
    end
end



function makeModelFrameCC(prs_path, pheno_path; header = 1)
    prs = CSV.read(prs_path; delim = "\t", header = header)

    fitted_path = joinpath(pheno_path, "results/casecontrol", basename(pheno_path), "fitted.auto.tsv")
    fitted = CSV.read(fitted_path; delim = "\t", header = [:PN, :fitted])

    case_path = joinpath(pheno_path, "results/casecontrol", basename(pheno_path), "aff.auto.pnlist")
    case = CSV.read(case_path; delim = "\t", header = [:PN])

    pns = Set(case.PN)

    m = join(prs, fitted, on = :PN)
    m[!,:affected] .= 0
    m.affected[in(pns).(m.PN)] .= 1

    m
end


"""
    estimate_h2(z, lds, n)

Estimate total heritability ``h^2`` for markers using LD score
regression with intercept forced to 1.

# Arguments

- `lds::Vector`: LD score as computed by `ldscore(LDMatrix)`
- `n::Integer`: effective number of samples (``4 n_0*n_1/(n_0+n_1)`` for case-control)
"""
function estimate_h2(z, lds, n)
    M = length(lds)
    meanlds = sum(lds) / M

    h2 = (mean(z .^2) - 1) / n / meanlds * M
end


"""
    estimate_neff(beta, sebeta, freq)

Estimate effective number of samples in case-control study.

Effective number of samples is total number of samples for a cohort
with same number of samples and controls, that would result in the
observed `sebeta`.

`freq` can be minor or major allele freq, 
as only `freq*(1-freq)` is being used.
"""
function estimate_neff(beta, sebeta, freq)
    neff = 2. ./ (sebeta .^2 .* freq .* (1.0 .- freq))
    median(neff[isfinite.(neff)])
end


"""
    snp_join(df1, df2; on = [:Chrom, :Pos],
             alleles1 = [:A1_1, :A2_1], alleles2 = [:A1_2, :A2_2],
             matchcol = :sign)

Join data frames `df1` and `df2` and match variants.

Data frames `df1` and `df2` are first matched by columns specified in
`on`. Then column `matchcol` is set to

- `matchcol == 1` if values for `alleles1` and `alleles2` match
- `matchcol == -1` if values for `alleles1` and `alleles2` are swapped
- `matchcol == 0` otherwise
"""
function snp_join(df1, df2; on = [:Chrom, :Pos],
                  alleles1 = [:A1_1, :A2_1], alleles2 = [:A1_2, :A2_2],
                  matchcol = :sign)

    m = innerjoin(df1, df2, on = on)
    m[!, matchcol] .= 0
    m[(m[!, alleles1[1]] .== m[!, alleles2[1]]) .& (m[!, alleles1[2]] .== m[!, alleles2[2]]), matchcol] .= 1
    m[(m[!, alleles1[1]] .== m[!, alleles2[2]]) .& (m[!, alleles1[2]] .== m[!, alleles2[1]]), matchcol] .= -1

    n_match = sum(m[!, matchcol] .== 1)
    n_mismatch = sum(m[!, matchcol] .== 0)
    n_flip = sum(m[!, matchcol] .== -1)

    println("Found $n_match matching, $n_flip flipped, and $n_mismatch mismatching SNPs")
    m
end


"""
read chi markers 

return data frame with columns
Chrom, Name, Pos, A1, A2
"""
function read_chi_markers(chidir::AbstractString, chrom)
    path = joinpath(chidir, chrom * ".markers")
    tmp = CSV.read(path; delim="\t")[!,1:5]
    rename!(tmp, [:Chrom_, :Name, :Pos, :A1, :A2])
    tmp[!,:Chrom] .= "chr" .* string.(tmp.Chrom_)
    tmp[!, [:Chrom, :Name, :Pos, :A1, :A2]]
end


function get_chi_markers(dir::AbstractString, markers::Vector)
    requested = in(Set(markers))
    files = sort([joinpath(dir, f) for f in readdir(dir) if endswith(f, ".markers")])
    tmp = CSV.read(files[1]; delim="\t")
    tmp = tmp[requested.(tmp.name), 1:5]
    rename!(tmp, [:Chrom, :Name, :Pos, :A1, :A2])
    res = tmp
    for f in files[2:end]
        tmp = CSV.read(f; delim="\t")
        tmp = tmp[requested.(tmp.name),1:5]
        rename!(tmp, [:Chrom, :Name, :Pos, :A1, :A2])
        res = vcat(res, tmp)
    end
    res
end

    
function postfix_names!(df::DataFrame, postfix::String)
    oldnames = [string(x) for x in names(df)]
    newnames = [ x * postfix for x in oldnames]
    rename!(df, newnames)
end
