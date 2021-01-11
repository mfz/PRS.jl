
#
# Implementation of PRS similar to lassosum
#
# Mak et al.(2018), Polygenic scores via penalized regression on
# summary statistics
#

export z2cor,
    elnetg!,
    elnetg_path



"""
    z2cor(z, n)

Convert GWAS z-value to phenotype-genotype correlation coefficient.
"""
z2cor(z, n) = z ./ sqrt.(n .- 1 .+ z.^2)



"""
    elnetg!(X, r, λ, s, β = zeros(Float64, length(r));
            maxiter = 10_000, thresh = 1e-4)

Solve elastic net via coordinate descent.

# Arguments

- `X::Matrix` nsubj x nmarkers normalized genotype matrix
              need to pass in ac*sqrt((1-s)/nsubj),
              where ac is column normalized (μ=0, σ=1) allele count matrix
- `r::Vector` nmarkers x 1 vector of correlation between phenotype and genotypes
- `λ::Real`   shrinkage parameter for 1-norm
- `s::Real`   shrinkage parameter for 2-norm (LD)
- `β::Vector` warm start (in) and result (out) for solution vector
- `maxiter::Int` maximum number of iterations
- `thresh::Real` maximum change in β to be called converged

"""
function elnetg!(X, r, λ, s, β;
                 maxiter = 10_000, thresh = 1e-4,
                 verbose = false)

    p = length(r) # number of markers
    @assert p == size(X, 2) "Size mismatch: length(r) = $p, size(X, 2) = $(size(X, 2))"
    
    sd = std(X; dims = 1)[:]
    diag = ones(Float64, p) * (1.0 - s)
    diag[sd .== 0.0] .= 0.0

    yhat = X * β

    converged = false
   
    for k = 1:maxiter
        Δ = 0.0
        for j = 1:p
            β_j = β[j]
            β[j] = 0.0
            u_j = diag[j]*β_j + r[j]  - yhat'*X[:,j]


            abs_u_j = abs(u_j)
            #println("j=$j, abs_u_j = $abs_u_j")
            if abs_u_j > λ
                β[j] = sign(u_j) * (abs_u_j - λ) / (diag[j] + s)
            end
            Δ_j = β[j] - β_j
            Δ = max(Δ, abs(Δ_j))
            yhat += Δ_j * X[:,j]
        end

        if Δ < thresh
            converged = true
            break
        end

        if verbose
            print("\rλ = $λ, s = $s: Δ = $Δ  ($k/$maxiter)")   
        end
        
    end
    obj = yhat'*yhat - 2*β'*r + 2*λ*sum(β) + s*β'*β
    
    if verbose
        println("\rλ = $λ, s = $s: obj = $obj, converged = $converged")
    end

    (converged = converged, obj = obj) 
end


"""
    elnetg_path(X, r, λs, s;
                maxiter = 10_000, thresh = 1e-4)

Solve elastic net for path along `λs` using warm starts.

# Arguments

- `X::Matrix` nsubj x nmarkers normalized genotype matrix
              column normalized (μ = 0, σ = 1)
- `r::Vector` nmarkers x 1 vector of correlation coefficients 
              between phenotype and genotypes
- `λs::Vector` shrinkage parameter path for 1-norm in decreasing order
- `s::Real`   shrinkage parameter for 2-norm (LD)
- `β::Vector` warm start (in) and result (out) for solution vector
- `maxiter::Int` maximum number of iterations
- `thresh::Real` maximum change in β to be called converged

"""
function elnetg_path(X, r, λs, s;
                     maxiter = 1000, thresh = 1e-4, verbose = false)

    nl = length(λs)
    p = length(r)

    β = zeros(Float64, p, nl)
    converged = zeros(Bool, nl)

    nsubj = size(X, 1)
    Xs = X .* sqrt((1.0 - s)/nsubj)
    ret = elnetg!(Xs, r, λs[1], s, view(β, :, 1);
                  maxiter = maxiter, thresh = thresh, verbose = verbose)
    converged[1] = ret.converged
    
    for l = 2:nl
        @assert λs[l] < λs[l-1] "λs not in decreasing order"
        β[:, l] .= β[:, l-1] # use warm start
        ret = elnetg!(Xs, r, λs[l], s, view(β, :, l);
                     maxiter = maxiter, thresh = thresh, verbose = verbose)
        converged[l] = ret.converged
        
    end
    
    (β = β,
     λs = λs,
     converged = converged,
     s = s)

end



#
# Implementations using sparse LD matrix instead of genotype matrix.
# These have convergence problems, though.
#

"""
    elnet(R, r, λ, s, β = zeros(Float64, length(r));
          maxiter = 10_000, thresh = 1e-4)

Solve elastic net via coordinate descent.

# Arguments

- `R::Matrix` p x p LD matrix
- `r::Vector` p x 1 vector of correlation between phenotype and genotypes
- `λ::Real`   shrinkage parameter for 1-norm
- `s::Real`   shrinkage parameter for 2-norm (LD)
- `β::Vector` warm start (in) and result (out) for solution vector
- `maxiter::Int` maximum number of iterations
- `thresh::Real` maximum change in β to be called converged

"""
function elnet!(R, r, λ, s, β;
                maxiter = 10_000, thresh = 1e-4,
                verbose = false)

    p = length(r)
    @assert p == size(β, 1) "Size mismatch: length(r) = $p, size(β, 1) = $(size(β, 1))"
    
    converged = false

    for k = 1:maxiter
        Δ = 0.0   # maximum change in β
        for j = 1:p
            β_j = β[j]
            β[j] = 0.0
            u_j = s ≈ 1 ? r[j] : (r[j] - (1.0-s) * β' * R[:,j]) 
            abs_u_j = abs(u_j)
            if abs_u_j > λ
                β[j] = sign(u_j) * (abs_u_j - λ) / (s + (1-s)*R[j,j])
            end
            Δ_j = β[j] - β_j
            Δ = max(Δ, abs(Δ_j))
        end

        if Δ < thresh
            converged = true
            break
        end

        if verbose
            print("\rλ = $λ, s = $s: Δ = $Δ  ($k/$maxiter)")   
        end

        
    end

    converged
end



"""
    elnet_path(R, r, λs, s;
               maxiter = 10_000, thresh = 1e-4)

Solve elastic net for path along λs using warm starts.

HAS CONVERGENCE PROBLEMS. USE elnetg_path INSTEAD

# Arguments

- `R::Matrix` p x p LD matrix
- `r::Vector` p x 1 vector of correlation coefficients 
              between phenotype and genotypes
- `λs::Vector` shrinkage parameter path for 1-norm in decreasing order
- `s::Real`   shrinkage parameter for 2-norm (LD)
- `β::Vector` warm start (in) and result (out) for solution vector
- `maxiter::Int` maximum number of iterations
- `thresh::Real` maximum change in β to be called converged

"""
function elnet_path(R, r, λs, s;
                    maxiter = 10_000, thresh = 1e-4)

    nl = length(λs)
    p = length(r)

    β = zeros(Float64, p, nl)
    converged = zeros(Bool, nl)
    
    converged[1] = elnet!(R, r, λs[1], s, view(β, :, 1);
                         maxiter = maxiter, thresh = thresh)
    
    for l = 2:nl
        @assert λs[l] < λs[l-1] "λs not in decreasing order"
        β[:, l] .= β[:, l-1] # use warm start
        converged[l] = elnet!(R, r, λs[l], s, view(β, :, l);
                              maxiter = maxiter, thresh = thresh)
    end
    
    (β = β,
     λs = λs,
     converged = converged,
     s = s)

end
