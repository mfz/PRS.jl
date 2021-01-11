# PRS.jl Documentation

```@contents
```

Implementation of ldpred, EB, and lassosum in Julia.

## PlinkReader

```@docs
PlinkReader(path; markerIndex = false, sampleIndex = false)

Sample
nsamples(p::PlinkReader)
samples(p::PlinkReader)
sample_index(p::PlinkReader, iid::String)

Marker
nmarkers(p::PlinkReader)
markers(p::PlinkReader)
marker_index(p::PlinkReader, id::String)
markersDF(p::PlinkReader)

getindex(p::PlinkReader, s::Int, m::Int)

dosageMatrix(p::PlinkReader, markerIdx, sampleIdx = nothing;
             normalize = true)
```

## LDMatrix

```@docs
LDMatrix(p::PlinkReader, mIdx0, mIdx1; alpha = 0.9, window = posWindow(1_000_000))
LDMatrix(path::AbstractString)

markers(ld::LDMatrix)
markersDF(ld::LDMatrix)

save(ld::LDMatrix, path::AbstractString)

ldscore(ld::LDMatrix)
```


## ldpred

```@docs
ldpred_gibbs(z, D, p, σ2, μ0; n_burnin = 100, n_iter = 500, verbose = false)

estimate_h2(z, lds, n)
estimate_neff(beta, sebeta, freq)
```


## lassosum

```@docs
z2cor(z, n)

elnetg!(X, r, λ, s, β = zeros(Float64, length(r));
        maxiter = 10_000, thresh = 1e-4)

elnetg_path(X, r, λs, s;
            maxiter = 10_000, thresh = 1e-4)
```


## utils

```@docs
snp_join(df1, df2;
         on = [], alleles1 = [], alleles2 = [], matchcol = :sign)
```
