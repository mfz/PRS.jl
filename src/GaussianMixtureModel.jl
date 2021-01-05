
export gmm,
    shrink,
    estimate_component_prob,
    estimate_var_μ,
    rand_gmm

using Distributions
using StatsBase

"""
fit Gaussian mixture model with all means set to 0.

A SNP contributing h2 to heritability has μ ~ N(0, var = N*h2)
and z|μ ~ N(μ, 1.0)

Here we decompose distribution of z into null component z ~ N(0; 1)
and components z ~ N(0; var = 1 + N*h2/M) where h2 is total heritability 

Due to population stratification etc. the null component might have
variance > 1. Then z needs to be rescaled such that variance of null component ≈ 1

dp0: Dirichlet pseudo counts for null component

"""
function gmm(x, π0, σ0; ϵ = 1e-3, maxiter = 1000, dp0 = 0, verbose = false, forcenull = false)

    n_data = length(x)
    n_cluster = length(π0)

    # Dirichlet prior
    dp = zeros(Float64, 1, n_cluster)
    dp[1] = dp0
    
    π = π0[:]'
    σ = σ0[:]'
    π_old = π0[:]'
    σ_old = σ0[:]'
    pc = zeros(Float64, n_data, n_cluster)

    llold = -Inf
    converged = false
    n_iter = 0
    
    for iter = 1:maxiter
        n_iter += 1
        # E-step
        # P(c_j|x_j) propto P(x_j|c_j) * P(c_j)
        for i = 1:n_data
            pc[i,:] = pdf.(Normal.(0., σ), x[i]) .* π  
        end # i
        pc .= pc ./ sum(pc; dims = 2)

        # M-step
        π .= (sum(pc; dims = 1) .+ dp) ./ (n_data + dp0)
        σ  .= sqrt.(sum(x.^2 .* pc; dims = 1) ./ sum(pc; dims = 1))
        if forcenull
            σ[1] = 1.0
        end
        
        # loglikelihood
        # sum_i log( sum_c P(x_i|c)*P(c) )
        ll = 0.0
        for i = 1:n_data
            ll += log(sum(pdf.(Normal.(0., σ), x[i]) .* π)) 
        end

        verbose && println("Iteration $iter : LogLikelihood $ll")
        
        if (maximum(abs.((π_old .- π) ./ π_old)) < ϵ) && (maximum(abs.((σ_old .- σ) ./ σ_old)) < ϵ)
            converged = true
            break
        end
        
        llold = ll
        π_old .= π
        σ_old = σ
        
    end
    
    (π = π, σ = σ, converged = converged, ll = llold, iter = n_iter)
end


"""
estimate component probability for z 

pc[i, k] : probability that z[i] belongs to component (π[k], σ[k])
"""
function estimate_component_prob(z, π, σ)
    @assert ndims(π) == ndims(σ) == 2 "π and σ need to be row vectors"
    @assert ndims(z) == 1 "z needs to be a column vector"
    @assert size(π) == size(σ) "π and σ need to have same dimensions"
    n_data = length(z)
    n_cluster = size(π, 2)
    pc = zeros(Float64, n_data, n_cluster)
    for i = 1:n_data
        pc[i,:] = pdf.(Normal.(0., σ), z[i]) .* π
    end
    pc ./ sum(pc; dims = 2)
end



function rand_gmm(π, σ)
    K = length(π)
    k = sample(1:K, Weights(π[:]))
    rand(Normal(0.0, σ[k]), 1)[1]
end



"""
compute shrinkage factor for z

i.e. estimate μ|z  where z|μ ~ N(μ; sd = 1)
and prior for μ given by  P(μ) = sum_k π[k] N(0; var = σ[k]^2 - 1)

σ2_μ for first component is set to 0
"""
function shrink(z, π, σ)

    pc = estimate_component_prob(z, π, σ)
    σ2_μ = max.(σ .^2 .- 1.0, 0.0)
    σ2_μ[1] = 0.0

    pc * (σ2_μ ./ (1.0 .+ σ2_μ))' 
end


"""
estimate var_μ = N h2 / M for each SNP

N  number of samples
h2 total heritability
M  number of SNPs

π: component probability as obtained from gmm
σ: standard deviation of component as obtained from gmm

1/var_μ needs to be added to diagonal of LD matrix 
"""
function estimate_var_μ(z, π, σ)
    
    pc = estimate_component_prob(z, π, σ)
    σ2_μ = max.(σ .^2 .- 1.0, 0.0)

    pc * σ2_μ'
end


function test_gmm()

    x = vcat(rand(Normal(0, 1.5), 100),
             rand(Normal(0,1.), 800),
             rand(Normal(0, 5.), 100))

    σ0 = [1., 2., 3.]
    π0 = [0.9, 0.05, 0.05]

    res = gmm(x, π0, σ0)
end
