var documenterSearchIndex = {"docs":
[{"location":"lassosum/#lassosum","page":"lassosum","title":"lassosum","text":"","category":"section"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"Mak et al.(2018), Polygenic scores via penalized regression on summary statistics","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"Consider linear regression (GWAS) problem","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"mathbfy = mathbfXbeta + mathbfepsilon","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"where mathbfX is the n x p genotype matrix (normalized to sd(cols) = 1) and mathbfy is n x 1 phenotype matrix (normalized to sd(cols) = 1)","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"Using regularization, we arrive at the LASSO with objective function","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"f(mathbfbeta) = (mathbfy - mathbfXbeta)^T (mathbfy - mathbfXbeta) + 2 lambda beta_1","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"f(mathbfbeta) = mathbfy^T mathbfy - 2 mathbfbeta^T mathbfX^T mathbfy + mathbfbeta^T mathbfX^T mathbfX mathbfbeta + 2 lambda beta_1","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"Observing that ","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"mathbfr = mathbfX^T mathbfy","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"is the genotype - phenotype correlation (obtained from p-value in summary statistics) and","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"mathbfR = mathbfX^T mathbfX","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"is the LD matrix, ","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"the objective can be rewritten as","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"f(mathbfbeta) = mathbfy^T mathbfy - 2 mathbfbeta^T mathbfr + mathbfbeta^T mathbfR mathbfbeta + 2 lambda beta_1","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"The LD matrix mathbfR could also be derived from reference data. In this case, mathbfr and mathbfR are not derived from the same genotype matrix mathbfX anymore, implying that the objective function is no longer a penalized least-squares (or LASSO) problem.","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"It can be turned into a LASSO problem when regularization","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"mathbfR_s = (1-s)mathbfR + smathbfI 0 le s le 1","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"is applied (see simple proof in paper). The shrinkage parameter s determines how much of the linkage disequilibrium between different variants is taken into account.","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"This renders the objective function an elastic net,","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"f(mathbfbeta) = mathbfy^T mathbfy - 2 mathbfbeta^T mathbfr + (1-s)mathbfbeta^T mathbfR mathbfbeta + s mathbfbeta^T mathbfbeta + 2 lambda beta_1","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":",","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"which can be solved using coordinate descent.","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"In each iteration of coordiante descent, the contribution of each single coordinate is optimized on its own. For a single coordinate beta_l, the contributions to f(beta) are","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"beta_l^2 (s + (1-s)mathbfR_ll) - 2 beta_l (r_l - (1-s)mathbfR_lbeta + (1-s)mathbfR_llbeta_l) + 2lambda beta_l","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"NOTE: Is not s + (1-s)mathbfR_ll = 1? Only if R_ll = 1, but for some markers the variance is 0, hence R_ll = 0","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"Defining ","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"u_l = r_l - (1-s)R_l beta + (1-s) R_ll beta_l","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"the minimum contribution is obtained at","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"beta_l = mathmrsgn(u_l) (u_l - lambda)  (s + (1-s)mathbfR_ll)  u_l - lambda  0","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"beta_l = 0 mathrmotherwise","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"NOTE: Is this not the same as is done in LDpred? I.e. we look at the contribution of all SNPs except l to mathbfr_l due to LD, and then apply some shrinkage to beta_l.","category":"page"},{"location":"lassosum/","page":"lassosum","title":"lassosum","text":"NOTE: Implementation where thresholded LD matrix is used suffers from divergence issues. Implementation using genotype matrix directly converges fine.","category":"page"},{"location":"#PRS.jl-Documentation","page":"Home","title":"PRS.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Implementation of ldpred, EB, and lassosum in Julia.","category":"page"},{"location":"#PlinkReader","page":"Home","title":"PlinkReader","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PlinkReader(path; markerIndex = false, sampleIndex = false)\n\nSample\nnsamples(p::PlinkReader)\nsamples(p::PlinkReader)\nsample_index(p::PlinkReader, iid::String)\n\nMarker\nnmarkers(p::PlinkReader)\nmarkers(p::PlinkReader)\nmarker_index(p::PlinkReader, id::String)\nmarkersDF(p::PlinkReader)\n\ngetindex(p::PlinkReader, s::Int, m::Int)\n\ngenotypeMatrix(p::PlinkReader, markerIdx, sampleIdx = nothing;\n               normalize = true)","category":"page"},{"location":"#PRS.PlinkReader-Tuple{Any}","page":"Home","title":"PRS.PlinkReader","text":"PlinkReader(path; markerIndex = false, sampleIndex = false)\n\nStructure to access PLINK .bed, .bim, and .fam files at path.\n\nArguments\n\nmarkerIndex::Bool should index marker -> idx be created?\nsampleIndex::Bool should index sample -> idx be created?\n\n\n\n\n\n","category":"method"},{"location":"#PRS.Sample","page":"Home","title":"PRS.Sample","text":"Sample(fid, iid, father, mother, sex)\n\nPLINK sample info\n\n\n\n\n\n","category":"type"},{"location":"#PRS.nsamples-Tuple{PlinkReader}","page":"Home","title":"PRS.nsamples","text":"nsamples(p::PlinkReader)\n\nNumber of samples in Plink file.\n\n\n\n\n\n","category":"method"},{"location":"#PRS.samples-Tuple{PlinkReader}","page":"Home","title":"PRS.samples","text":"samples(p::PlinkReader)\n\nReturn Array of Sample structs for samples in Plink file.\n\n\n\n\n\n","category":"method"},{"location":"#PRS.sample_index-Tuple{PlinkReader,String}","page":"Home","title":"PRS.sample_index","text":"sample_index(p::PlinkReader, iid::String)\n\nReturn index of sample iid in PlinkReader. Uses index if available, otherwise linear search. \n\n\n\n\n\n","category":"method"},{"location":"#PRS.Marker","page":"Home","title":"PRS.Marker","text":"Marker(chrom, id, cm, pos, a1, a2)\n\nPLINK marker info\n\n\n\n\n\n","category":"type"},{"location":"#PRS.nmarkers-Tuple{PlinkReader}","page":"Home","title":"PRS.nmarkers","text":"nmarkers(p::PlinkReader)\n\nNumber of markers in Plink file.\n\n\n\n\n\n","category":"method"},{"location":"#PRS.markers-Tuple{PlinkReader}","page":"Home","title":"PRS.markers","text":"markers(p::PlinkReader)\n\nReturn array of Marker structs for markers in Plink file.\n\n\n\n\n\n","category":"method"},{"location":"#PRS.marker_index-Tuple{PlinkReader,String}","page":"Home","title":"PRS.marker_index","text":"marker_index(p::PlinkReader, id::String)\n\nReturn index of marker id in PlinkReader. Uses index if available, otherwise linear search.\n\n\n\n\n\n","category":"method"},{"location":"#PRS.markersDF-Tuple{PlinkReader}","page":"Home","title":"PRS.markersDF","text":"markersDF(p::PlinkReader)\n\nGet data frame with marker info.\n\nReturns data frame with columns\n\nChrom\nName\ncM\nPos\nA1\nA2\nIdx\n\n\n\n\n\n","category":"method"},{"location":"#Base.getindex-Tuple{PlinkReader,Int64,Int64}","page":"Home","title":"Base.getindex","text":"getindex(p::PlinkReader, s::Int, m::Int)\n\nRetrieve genotype of sample s and marker m.\n\nPlink convention for representation of genotypes is as follows:\n\n0b00 Hom1\n0b01 Het \n0b10 missing\n0b11 Hom2\n\n\n\n\n\n","category":"method"},{"location":"#PRS.genotypeMatrix","page":"Home","title":"PRS.genotypeMatrix","text":"genotypeMatrix(p::PlinkReader, markerIdx, sampleIdx = nothing;\n               normalize = true)\n\nCreate genotype matrix (samples x markers) for markerIdx and sampleIdx from PlinkReader. Here genotype is expressed as count of alternative allele (not as in Plink convention).\n\nReplace missing genotypes with mean genotype. If normalize = true, normalize markers to mean μ = 0 and standard deviation σ = 1.\n\n\n\n\n\n","category":"function"},{"location":"#LDMatrix","page":"Home","title":"LDMatrix","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"LDMatrix(p::PlinkReader, mIdx0, mIdx1; alpha = 0.9, window = posWindow(1_000_000))\nLDMatrix(path::AbstractString)\n\nmarkers(ld::LDMatrix)\nmarkersDF(ld::LDMatrix)\n\nsave(ld::LDMatrix, path::AbstractString)\n\nldscore(ld::LDMatrix)","category":"page"},{"location":"#PRS.LDMatrix-Tuple{PlinkReader,Any,Any}","page":"Home","title":"PRS.LDMatrix","text":"LDMatrix(p::PlinkReader, mIdx0, mIdx1; alpha = 0.9, window = posWindow(1_000_000))\n\nCompute LDMatrix from genotypes in Plink file.\n\n\n\n\n\n","category":"method"},{"location":"#PRS.LDMatrix-Tuple{AbstractString}","page":"Home","title":"PRS.LDMatrix","text":"LDMatrix(path::AbstractString)\n\nLoad LDMatrix from .bim and .lds files at path.\n\n\n\n\n\n","category":"method"},{"location":"#PRS.markers-Tuple{LDMatrix}","page":"Home","title":"PRS.markers","text":"markers(ld::LDMatrix)\n\nGet array of Marker structs for all markers in LDMatrix ld.\n\n\n\n\n\n","category":"method"},{"location":"#PRS.markersDF-Tuple{LDMatrix}","page":"Home","title":"PRS.markersDF","text":"markersDF(ld::LDMatrix)\n\nGet data frame with marker info.\n\nReturns data frame with columns\n\nChrom\nName\ncM\nPos\nA1\nA2\nIdx\n\n\n\n\n\n","category":"method"},{"location":"#PRS.save-Tuple{LDMatrix,AbstractString}","page":"Home","title":"PRS.save","text":"save(ld::LDMatrix, path::AbstractString)\n\nSave LDMatrix to path.\n\n\n\n\n\n","category":"method"},{"location":"#PRS.ldscore-Tuple{LDMatrix}","page":"Home","title":"PRS.ldscore","text":"ldscore(ld::LDMatrix)\n\nCompute LD scores for markers in LDMatrix ld.\n\n\n\n\n\n","category":"method"},{"location":"#ldpred","page":"Home","title":"ldpred","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ldpred_gibbs(z, D, p, σ2, μ0; n_burnin = 100, n_iter = 500, verbose = false)\n\nestimate_h2(z, lds, n)\nestimate_neff(beta, sebeta, freq)","category":"page"},{"location":"#PRS.ldpred_gibbs-NTuple{5,Any}","page":"Home","title":"PRS.ldpred_gibbs","text":"ldpred_gibbs(z, D, p, σ2, μ0; n_burnin = 100, n_iter = 500, verbose = false)\n\nRun LDpred Gibbs sampler.\n\nArguments\n\nz::Vector: z-scores\nD::LDMatrix: LD matrix\np::Real : proportion of variants deemed to be causal\nσ2::Real: Nh²M as estimated from LDscore regression (estimate_h2, or mean(z^2)/mean(lds))             Note: prior variance of non-null component of μ, σ2_μ = σ2p \nμ0::Vector: starting estimate (e.g. from infinite model)\n\n\n\n\n\n","category":"method"},{"location":"#PRS.estimate_h2-Tuple{Any,Any,Any}","page":"Home","title":"PRS.estimate_h2","text":"estimate_h2(z, lds, n)\n\nEstimate total heritability h^2 for markers using LD score regression with intercept forced to 1.\n\nArguments\n\nlds::Vector: LD score as computed by ldscore(LDMatrix)\nn::Integer: effective number of samples (4 n_0*n_1(n_0+n_1) for case-control)\n\n\n\n\n\n","category":"method"},{"location":"#PRS.estimate_neff-Tuple{Any,Any,Any}","page":"Home","title":"PRS.estimate_neff","text":"estimate_neff(beta, sebeta, freq)\n\nEstimate effective number of samples in case-control study.\n\nEffective number of samples is total number of samples for a cohort with same number of samples and controls, that would result in the observed sebeta.\n\nfreq can be minor or major allele freq,  as only freq*(1-freq) is being used.\n\n\n\n\n\n","category":"method"},{"location":"#lassosum","page":"Home","title":"lassosum","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"z2cor(z, n)\n\nelnetg!(X, r, λ, s, β = zeros(Float64, length(r));\n        maxiter = 10_000, thresh = 1e-4)\n\nelnetg_path(X, r, λs, s;\n            maxiter = 10_000, thresh = 1e-4)","category":"page"},{"location":"#PRS.z2cor-Tuple{Any,Any}","page":"Home","title":"PRS.z2cor","text":"z2cor(z, n)\n\nConvert GWAS z-value to phenotype-genotype correlation coefficient.\n\n\n\n\n\n","category":"method"},{"location":"#PRS.elnetg_path-NTuple{4,Any}","page":"Home","title":"PRS.elnetg_path","text":"elnetg_path(X, r, λs, s;\n            maxiter = 10_000, thresh = 1e-4)\n\nSolve elastic net for path along λs using warm starts.\n\nArguments\n\nX::Matrix nsubj x nmarkers normalized genotype matrix             column normalized (μ = 0, σ = 1)\nr::Vector nmarkers x 1 vector of correlation coefficients              between phenotype and genotypes\nλs::Vector shrinkage parameter path for 1-norm in decreasing order\ns::Real   shrinkage parameter for 2-norm (LD)\nβ::Vector warm start (in) and result (out) for solution vector\nmaxiter::Int maximum number of iterations\nthresh::Real maximum change in β to be called converged\n\n\n\n\n\n","category":"method"},{"location":"#utils","page":"Home","title":"utils","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"snp_join(df1, df2;\n         on = [], alleles1 = [], alleles2 = [], matchcol = :sign)","category":"page"},{"location":"#PRS.snp_join-Tuple{Any,Any}","page":"Home","title":"PRS.snp_join","text":"snp_join(df1, df2; on = [:Chrom, :Pos],\n         alleles1 = [:A1_1, :A2_1], alleles2 = [:A1_2, :A2_2],\n         matchcol = :sign)\n\nJoin data frames df1 and df2 and match variants.\n\nData frames df1 and df2 are first matched by columns specified in on. Then column matchcol is set to\n\nmatchcol == 1 if values for alleles1 and alleles2 match\nmatchcol == -1 if values for alleles1 and alleles2 are swapped\nmatchcol == 0 otherwise\n\n\n\n\n\n","category":"method"}]
}
