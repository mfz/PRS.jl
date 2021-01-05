

export LDMatrix,
    save,
    posWindow,
    cMWindow,
    markers,
    markersDF,
    ldscore

using Distributions
using LinearAlgebra
using SparseArrays

cMWindow(cM) = (m1, m2) -> abs(m1.cm - m2.cm) <= cM
posWindow(dist) = (m1, m2) -> abs(m1.pos - m2.pos) <= dist


struct LDMatrix
    markers::Vector{Marker}
    markerIndex::Dict{String, Int64}
    r::Symmetric{Float64, SparseMatrixCSC{Float64, Int64}}
end


"""
    LDMatrix(p::PlinkReader, mIdx0, mIdx1; alpha = 0.9, window = posWindow(1_000_000))

Compute LDMatrix from genotypes in Plink file.
"""
function LDMatrix(p::PlinkReader, mIdx0, mIdx1; alpha = 0.9, window = posWindow(1_000_000))
    LDMatrix(markers(p)[mIdx0:mIdx1],
             Dict([m.id => idx for (idx, m) in enumerate(markers(p))]),
             Symmetric(ld(p, mIdx0, mIdx1; alpha = alpha, window = window), :L))
end

"""
    LDMatrix(path::AbstractString)

Load LDMatrix from .bim and .lds files at `path`.
"""
function LDMatrix(path::AbstractString)
    markers = load_bim(path)
    r = loadSparseCSC(path * ".lds")
    LDMatrix(markers,
             Dict([m.id => idx for (idx, m) in enumerate(markers)]),
             Symmetric(r, :L))
end

"""
    markers(ld::LDMatrix)

Get array of `Marker` structs for all markers in LDMatrix `ld`.
"""
markers(ld::LDMatrix) = ld.markers
marker_index(ld::LDMatrix, marker) = ld.markerIndex[marker]
Base.size(ld::LDMatrix) = size(ld.r)


"""
    markersDF(ld::LDMatrix)

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
function markersDF(ld::LDMatrix)
    m = DataFrame(Chrom = String[], Name = String[], cM = Float64[],
                  Pos = Int64[], A1 = String[], A2 = String[],
                  Idx = Int64[])

    for (idx, mrk) in enumerate(markers(ld))
        push!(m, (mrk.chrom, mrk.id, mrk.cm, mrk.pos, mrk.a1, mrk.a2, idx))
    end
    m
end


"""
    corr_thresh(n, alpha = 0.05)

Compute significance threshold for Pearson correlation.

r (estimated from n data pairs) significant at alpha
if r >= t/sqrt(n-2+t*t)
where t is quantile(TDist(n-2), 1-alpha/2)

see e.g. Wikipedia article for Pearson correlation
"""
function corr_thresh(n, alpha = 0.05)
    n > 2 || return Inf
    t = quantile(TDist(n - 2), 1.0 - alpha/2.0)
    t / sqrt(n - 2 + t*t)
end



# 20201222 set diagonal element to 0.0 if variance = 0.0


"""
    ld(p::PlinkReader, mIdx0, mIdx1; alpha = 0.9, window = posWindow(1_000_000))

Compute LD block for markers mIdx0..mIdx1 (usually one chromosome), 
restricting to window and r significant at alpha.

output is sparse LD matrix with indices 1 .. (mIdx1 - mIdx0 + 1)
"""
function ld(p::PlinkReader, mIdx0, mIdx1; alpha = 0.9, window = posWindow(1_000_000))

    # a2counts[gt + 1] = counts of allele2 for gt
    # Plink convention is
    # 0b00 Hom1, 0b01 Het, 0b10 missing, 0b11 Hom2
    a2counts = [0, 1, -1, 2]

    thresh = corr_thresh.(1:nsamples(p), alpha)
    

    colInd = Int64[1]    # colInd[j] index into rowInd, where column j starts
    rowInd = Int64[]
    nzvalue = Float64[]  # make this Float64 so we can use LinAlg routines from SparseSuite  

    
    for m0 = mIdx0:mIdx1
        # compute variance of marker m0
        xSum0 = 0
        xxSum0 = 0
        nGood = 0
        for s = 1:nsamples(p)
            x = a2counts[p[s, m0] + 1]
            x < 0 && continue
            nGood += 1
            xSum0 += x
            xxSum0 += x*x
        end # s
        #println("m0 = $m0, xSum0 = $xSum0, xxSum0 = $xxSum0, nGood = $nGood, var=$(xxSum0 - xSum0*xSum0/nGood)")
        # set diagonal to 1f0
        # ONLY DO THIS IF VARIANCE > 0!
        push!(rowInd, m0 - mIdx0 + 1)
        if (xxSum0 - xSum0*xSum0/nGood > 0.0) 
            push!(nzvalue, 1f0)
               
            for m =(m0+1):mIdx1
                nGood = 0  # number of valid x and y genotypes
                xSum = xSum0
                xxSum = xxSum0
                ySum = 0
                yySum = 0
                xySum = 0
                
                # check we are within window, or break
                window(markers(p)[m], markers(p)[m0])||break
                
                
                for s = 1:nsamples(p)
                    x = a2counts[p[s, m0] + 1]
                    x == -1 && continue
                    
                    y = a2counts[p[s, m] + 1]
                    if y == -1
                        # here y genotype missing but x genotype present
                        # need to correct xSum and xxSum for missing pair
                        xSum -= x
                        xxSum -= x*x
                        continue
                    end

                    nGood += 1
                    ySum += y
                    yySum += y*y
                    xySum += x*y
                    
                end # s

                nGood == 0 && continue
                
                cov_xy = xySum - xSum*ySum/nGood
                var_x = xxSum - xSum*xSum/nGood
                var_y = yySum - ySum*ySum/nGood

                r = cov_xy / sqrt(var_x * var_y)
                
                if abs(r) >= thresh[nGood]
                    push!(rowInd, m - mIdx0 + 1)
                    push!(nzvalue, r)
                end
                
            end # m
        
        else
            push!(nzvalue, 0f0)
        end #Var(m0) > 0
        push!(colInd, length(rowInd) + 1)
    end # m0
    
    # create sparse matrix from colInd, rowInd, value
    SparseMatrixCSC(mIdx1 - mIdx0 + 1,
                    mIdx1 - mIdx0 + 1,
                    colInd,
                    rowInd,
                    nzvalue)

end


function writeSparseCSC( s::SparseMatrixCSC, path::AbstractString)

    fh = open(path, "w")
    
    nb = 0
    
    m = s.m
    n = s.n
    nnz = length(s.nzval)

    nb += write(fh, Int64(m), Int64(n), Int64(nnz))

    for i = 1:(n+1)
        nb += write(fh, Int64(s.colptr[i]))
    end
    
    for i = 1:nnz
        nb += write(fh, Int64(s.rowval[i]))
    end

    for i = 1:nnz
        nb += write(fh, Float64(s.nzval[i]))
    end

    close(fh)
    
end

function loadSparseCSC(path::AbstractString)

    fh = open(path, "r")
    
    m = read(fh, Int64)
    n = read(fh, Int64)
    nnz = read(fh, Int64)

    offset = position(fh)
    colptr = Mmap.mmap(fh, Vector{Int64}, n + 1, offset; shared = false)

    offset += (n + 1) * 8
    rowval = Mmap.mmap(fh, Vector{Int64}, nnz, offset; shared = false)

    offset += nnz * 8
    nzval = Mmap.mmap(fh, Vector{Float64}, nnz, offset; shared = false)

    SparseMatrixCSC(m, n, colptr, rowval, nzval)
    
end


"""
    save(ld::LDMatrix, path::AbstractString)

Save `LDMatrix` to `path`.
"""
function save(ld::LDMatrix, path::AbstractString)
    save_bim(ld.markers, path * ".bim")
    writeSparseCSC(ld.r.data, path * ".lds")
end


function save_bim(markers, path)
    fh = open(path, "w")
    for m in markers
        println(fh, "$(m.chrom)\t$(m.id)\t$(m.cm)\t$(m.pos)\t$(m.a1)\t$(m.a2)")
    end
    close(fh)
end

"""
    ldscore(ld::LDMatrix)

Compute LD scores for markers in LDMatrix `ld`.
"""
function ldscore(ld::LDMatrix)
    # get underlying SparseMatrixCSC
    s = ld.r.data
    ldsc = diag(s)   #ones(Float64, s.n)  # diagonal elts can be 0.0
    
    for col = 1:(s.n-1)
        for idx = s.colptr[col]:(s.colptr[col+1]-1)
            row = s.rowval[idx]
            x = s.nzval[idx]
            # x is now r[row, col]
            row != col || continue
            r2 = x*x
            ldsc[col] += r2
            ldsc[row] += r2
        end
       
    end
    ldsc
end
