# Tutorial

Download data from "Basic Tutorial for Polygenic Risk Score Analyses" at https://choishingwan.github.io/PRS-Tutorial.

## lassosum

Load required packages and set paths to data files.

```julia
using DataFrames
using CSV
using Distributions: Normal, cquantile
using PRS

DATA_DIR = "/Users/florian/workAtHome/PRSdata"
SSF_PATH = joinpath(DATA_DIR, "Height.QC")
PLINK_PATH = joinpath(DATA_DIR, "EUR.QC")

RESULT_DIR = mkdir("/tmp/test_PRS_lassosum")
```

Load the summary data in a data frame and add "chr" prefix to chromosomes.

```julia
ssf = DataFrame!(CSV.File(SSF_PATH; delim="\t"))
ssf.Chrom = "chr" .* string.(ssf.CHR)
ssf.Pos = ssf.BP
```

The lassosum algorithm works on genotype-phenotype
correlations. Convert p-value to those.

```julia
ssf.Z = sign.(log.(ssf.OR)) .* cquantile.(Normal(),ssf.P / 2.0)
ssf.COR = z2cor.(ssf.Z, ssf.N)
```

```julia
plink = PlinkReader(PLINK_PATH)
m = markersDF(plink)
postfix_names!(m, "_GT")
```

Remove markers with mismatching alleles. 
TODO: remove ambiguous alleles as well

```julia
data = snp_join(ssf, m, on = [:Chrom => :Chrom_GT, :BP => :Pos_GT],
	       alleles1 = [:A1, :A2],
	       alleles2 = [:A1_GT, :A2_GT],
	       matchcol = :sign_GT)

data = data[data.sign_GT .!= 0, :]
```

Run lassosum on chromosome blocks. Note that the original lassosum
algorithm is run on LD blocks (or uniform blocks) instead. Need to
assess performance difference here.


```julia
ss = [0.2, 0.5, 0.9, 1.]
λs = exp.(range(log(0.01), log(0.001); length = 10)) # sorted from large to small

chrom = "chr21"

println("Processing chromosome $chrom")

cols = [:Chrom, :Pos, :A1, :A2]
block = data[data.Chrom .== chrom,:]
r = block.COR .* block.sign_GT
X = dosageMatrix(plink, block.Idx_GT)
for sidx in 1:length(ss)
  println("s = $(ss[sidx])")
  res = elnetg_path(X, r, λs, ss[sidx]; verbose=true)

  for lidx in 1:length(λs)
    col = Symbol("eff_lassosum_$(sidx)_$(lidx)")
    block[!, col] .= res.β[:,lidx] .* block.sign_GT
    push!(cols, col)
  end
  CSV.write(joinpath(RESULT_DIR, "$chrom.weights"), b[!, cols]; delim = "\t")
end	
```
