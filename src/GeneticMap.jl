export GeneticMap


using Interpolations

const DECODE_GENETIC_MAP = "/odinn/data/extdata/deCODE/GenMap/20150202_hg38/decode.genmap_hg38.gor"

abstract type AbstractGeneticMap end

struct GeneticMap <: AbstractGeneticMap
    path::String
    interp::Dict
end


function GeneticMap(path = DECODE_GENETIC_MAP)
    chromDict = Dict()
    fh = open(path)
    for line in eachline(fh)
        chrom, pos, name, cM = split(chomp(line), "\t")

        if !haskey(chromDict, chrom)
            chromDict[chrom] = (Float64[], Float64[])
        end

        push!(chromDict[chrom][1], parse(Float64, pos))
        push!(chromDict[chrom][2], parse(Float64, cM))
    end

    interp = Dict()
    for (chrom, pos_cM) in chromDict
        interp[chrom] = LinearInterpolation(pos_cM[1], pos_cM[2]; extrapolation_bc = Flat())
    end
    close(fh)
    
    GeneticMap(path, interp)
end


Base.show(io::IO, gm::GeneticMap) = print(io, "(GeneticMap: ", gm.path, ")")

(gm::GeneticMap)(chrom, pos) = gm.interp[chrom](pos)




