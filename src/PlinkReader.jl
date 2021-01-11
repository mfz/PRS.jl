
using Mmap
using DataFrames
using Statistics: mean, std

export Marker,
    Sample,
    PlinkReader,
    samples,
    markers,
    markersDF,
    nsamples,
    nmarkers,
    sample_index,
    marker_index,
    getindex,
    applyGeneticMap!,
    dosageMatrix

"""
    Marker(chrom, id, cm, pos, a1, a2)

PLINK marker info
"""
struct Marker
    chrom::String
    id::String
    cm::Float64
    pos::Int64
    a1::String
    a2::String
end


"""
    load_bim(path)

Load marker info from PLINK .bim file, adding 'chr' prefix if missing.
Return `Vector{Marker}`.
"""
function load_bim(path)
    
    markers = Marker[]
    bimpath = endswith(path, ".bim") ? path : path * ".bim"
    fh = open(bimpath, "r")
    for line in eachline(fh)
        cols = split(chomp(line))
        push!(markers, Marker(startswith(cols[1], "chr") ? cols[1] : "chr"*cols[1],
                              cols[2],
                              parse(Float64, cols[3]),
                              parse(Int64, cols[4]),
                              cols[5],
                              cols[6]))
    end
    close(fh)

    markers
end


"""
    Sample(fid, iid, father, mother, sex)

PLINK sample info
"""
struct Sample
    fid::String
    iid::String
    father::String
    mother::String
    sex::UInt8
end


"""
    load_fam(path)

Load sample info from PLINK .fam file. Return `Vector{Sample}`.
"""
function load_fam(path)

    samples = Sample[]
    fampath = endswith(path, ".fam") ? path : path * ".fam"
    fh = open(fampath, "r")
    for line in eachline(fh)
        cols = split(chomp(line))
        push!(samples, Sample(cols[1],
                              cols[2],
                              cols[3],
                              cols[4],
                              parse(UInt8, cols[5])))
    end
    close(fh)

    samples
end


"""
PLINK files (.bim, .fam, .bed)
"""
mutable struct PlinkReader
    path::String
    nmarkers::Int64
    nsamples::Int64
    markers::Vector{Marker}
    markerIndex::Dict{String, Int64}
    samples::Vector{Sample}
    sampleIndex::Dict{String, Int64}
    data::Array{UInt8, 2}
end


"""
    PlinkReader(path; markerIndex = false, sampleIndex = false)

Structure to access PLINK .bed, .bim, and .fam files at `path`.

# Arguments

- `markerIndex::Bool` should index marker -> idx be created?
- `sampleIndex::Bool` should index sample -> idx be created?
"""
function PlinkReader(path; markerIndex = false, sampleIndex = false)

    samples = load_fam(path)
    markers = load_bim(path)
    
    nsamples = length(samples)
    nmarkers = length(markers)

    fh = open(path * ".bed", "r")
    bedheader = read(fh, 3)

    @assert bedheader[1] == 0b01101100 && bedheader[2] == 0b00011011 && bedheader[3] == 0b01 "Only snp major BED v1.0 supported"

    #plinkdata = Mmap.mmap(fh, BitArray, (2, 4*ceil(Int64, 0.25 * nsamples), nmarkers))
    # see source code for Mmap.mmap(io, BitArray, dims)
    nrows = (nsamples + 3) >> 2
    data = Mmap.mmap(fh, Vector{UInt8}, nrows * nmarkers)
    
    PlinkReader(path,
                nmarkers,
                nsamples,
                markers,
                markerIndex ? Dict([m.id => idx for (idx, m) in enumerate(markers)]) : Dict(),
                samples,
                sampleIndex ? Dict([s.iid => idx for (idx, s) in enumerate(samples)]) : Dict(),
                reshape(data, (nrows, nmarkers)))
end


function Base.show(io::IO, p::PlinkReader)
    print(io, "<PLINK file ($(p.nsamples) samples x $(p.nmarkers) markers) at $(p.path)>")
end


"""
    samples(p::PlinkReader)

Return Array of `Sample` structs for samples in Plink file.
"""
samples(p::PlinkReader) = p.samples

"""
    markers(p::PlinkReader)

Return array of `Marker` structs for markers in Plink file.
"""
markers(p::PlinkReader) = p.markers

"""
    nsamples(p::PlinkReader)

Number of samples in Plink file.
"""
nsamples(p::PlinkReader) = p.nsamples

"""
    nmarkers(p::PlinkReader)

Number of markers in Plink file.
"""
nmarkers(p::PlinkReader) = p.nmarkers


"""
    markersDF(p::PlinkReader)

Get data frame with marker info.

Returns data frame with columns

- Chrom
- Name
- cM
- Pos
- A1
- A2
- Idx
"""
function markersDF(p::PlinkReader)
    m = DataFrame(Chrom = String[], Name = String[], cM = Float64[],
                  Pos = Int64[], A1 = String[], A2 = String[],
                  Idx = Int64[])

    for (idx, mrk) in enumerate(markers(p))
        push!(m, (mrk.chrom, mrk.id, mrk.cm, mrk.pos, mrk.a1, mrk.a2, idx))
    end
    m
end



"""
    sample_index(p::PlinkReader, iid::String)

Return index of sample `iid` in `PlinkReader`.
Uses index if available, otherwise linear search. 
"""
function sample_index(p::PlinkReader, iid::String)

    if length(p.sampleIndex) > 0
        return p.sampleIndex[iid]
    end
    
    s = samples(p)
    for i = 1:length(s)
        if s[i].iid == iid
            return i
        end
    end
    error("No such sample: $iid")
end


"""
    marker_index(p::PlinkReader, id::String)

Return index of marker `id` in `PlinkReader`.
Uses index if available, otherwise linear search.
"""
function marker_index(p::PlinkReader, id::String)

    if length(p.markerIndex) > 0
        return p.markerIndex[id]
    end
    
    m = markers(p)
    for i = 1:length(m)
        if m[i].id == id
            return i
        end
    end
    error("No such marker: $iid")
end    

"""
    getindex(p::PlinkReader, s::Int, m::Int)

Retrieve genotype of sample `s` and marker `m`.

Plink convention for representation of genotypes is as follows:

- `0b00` Hom1
- `0b01` Het 
- `0b10` missing
- `0b11` Hom2
"""
@inline function Base.getindex(p::PlinkReader, s::Int, m::Int)
    sp3 = s + 3
    (p.data[sp3 >> 2, m] >> ((sp3 & 0x03) << 1)) & 0x03
end

Base.size(p::PlinkReader) = (p.nsamples, p.nmarkers)


abstract type AbstractGeneticMap end

function applyGeneticMap!(gm::AbstractGeneticMap, p::PlinkReader)

    function applyGeneticMap(gm::GeneticMap, m::Marker)
        Marker(m.chrom,
               m.id,
               gm(m.chrom, m.pos),
               m.pos,
               m.a1,
               m.a2)
    end

    p.markers = [applyGeneticMap(gm, m) for m in p.markers]
    p
end


"""
    dosageMatrix(p::PlinkReader, markerIdx, sampleIdx = nothing;
                   normalize = true)

Create dosage matrix (samples x markers) for `markerIdx` and
`sampleIdx` from `PlinkReader`. Here dosage is expressed as count of
alternative allele.

Replace missing genotypes with mean dosage. If `normalize = true`,
normalize markers to mean μ = 0 and standard deviation σ = 1.
"""
function dosageMatrix(p::PlinkReader, markerIdx, sampleIdx = nothing;
                        normalize = true)

    sampleIdx == nothing && (sampleIdx = 1:nsamples(p))
    
    nm = length(markerIdx)
    ns = length(sampleIdx)
    g = zeros(Float64, ns, nm)

    # a2counts[gt + 1] = counts of allele2 for gt
    # Plink convention is
    # 0b00 Hom1, 0b01 Het, 0b10 missing, 0b11 Hom2
    a2counts = [0, 1, -1, 2]

    for (midx, m) in enumerate(markerIdx)
        data = view(p.data,:, m)
        counts = 0
        total = 0
        miss = 0
        for (sidx, s) in enumerate(sampleIdx)
            sp3 = s + 3
            plgt = data[sp3 >> 2] >> ((sp3 & 0x03) << 1) & 0x03

            a2count = a2counts[plgt + 1]
            g[sidx, midx] = a2count
            
            if a2count == -1
                miss += 1
            else
                counts += a2count
                total += 1
            end
        end

        # replace missing genotype with mean
        if miss > 0
            freq = counts/float(total)
            for s = 1:ns
                if g[s, midx] < 0.0
                    g[s, midx] = freq
                end
            end
        end

    end

    if normalize
        sd = std(g; dims = 1)
        sd[sd .== 0.0] .= 1.0 # avoid division by zero
        g ./= sd
        g .-= mean(g; dims = 1)
    end

    g

end
