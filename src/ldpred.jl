#
# Implementation of LDpred Gibbs sampler
#

#
# true effects μ distributed as
#
#        N(0, Nh^2/Mp)   with probability p
# μ_i ~  
#        0 with probability 1-p
#
# define σ2 := Nh^2/M
# Nh²/M can be obtained from LD-score regression
#
# σ2_μ = σ2/p

#
# observed (marginal) z-scores distributed as
#
# z_i|μ_i = N(μ_i, 1)
#
# σ2_z = 1 + σ2_μ
# σ_z = sqrt(1 + σ2_μ)

#
# residualized marginal effect
# i.e. correct SNP for effect of all other SNPs in LD
#
# z̃ = ẑ - (Dμ - μ)
#

#
# probability that variant j is causal
# is probability that z̃_i ~ N(0, 1 + Nh^2/Mp) vs z̃_i ~  N(0, 1)
#
#
# p_j = 1 / (1 + (1-p)/p * sqrt(1+σ2_μ) * exp(-0.5 * z̃^2 * σ2_μ / (1 + σ2_μ)))
#

#
# Posterior for μ_j is then
#
#             N(σ2_μ/(1 + σ2_μ) z̃_j, σ2_μ / (1 + σ2_μ))  with probability p_j  
# μ_j | z̃_j ~
#             0                                      with probability 1 - p_j
#

#
# posterior mean of μ_j is then
#
# ω_j = p_j * z̃_j * σ^2 / (1 + σ^2)
#


#
# Gibbs sampler
#
# for iter = 1:(n_burnin + n_iter)
#   for j = 1:n_variants
#      compute z̃_j
#      compute p_j
#      sample μ_j
#      compute ω_j
#      if iter > n_burnin
#        Ω_j += ω_j
#      end
#   end
# end
# μ = Ω / n_iter
#

using StatsBase

export ldpred_gibbs,
    ldpred_gibbs_EB,
    ldpred_gibbs_acc

"""
    ldpred_gibbs(z, D, p, σ2, μ0; n_burnin = 100, n_iter = 500, verbose = false)

Run LDpred Gibbs sampler.

# Arguments

- `z::Vector`: z-scores
- `D::LDMatrix`: LD matrix
- `p::Real` : proportion of variants deemed to be causal
- `σ2::Real`: ``Nh²/M`` as estimated from LDscore regression (`estimate_h2`, or `mean(z^2)/mean(lds)`)
              Note: prior variance of non-null component of ``μ``, ``σ2_μ = σ2/p`` 
- `μ0::Vector`: starting estimate (e.g. from infinite model)
"""
function ldpred_gibbs(z, D, p, σ2, μ0; n_burnin = 100, n_iter = 500, verbose = false)

    n_variants = length(z)
    μ = μ0[:]
    Ω = zeros(Float64, n_variants)

    σ2_μ = σ2 / p
    σ_z = sqrt(1.0 + σ2_μ)
    η = σ2_μ / (1.0 + σ2_μ)  # shrinkage factor
    
    for iter = 1:(n_burnin + n_iter)
        verbose && println("Iteration $iter/$(n_burnin + n_iter)")
        for j = 1:n_variants
            z̃_j = z[j] - μ'*D[:,j] + μ[j]
            p_j = 1.0 / (1.0 + (1.0-p)/p * σ_z * exp(-0.5 * z̃_j * z̃_j * η))
            μ[j] = (rand() <= p_j) ?  rand(Normal(η * z̃_j, sqrt(η))) : 0.0
            if iter > n_burnin
                Ω[j] += p_j * z̃_j * η
            end
        end
    end

    res = Ω ./ n_iter
    res[(!isfinite).(res)] .= 0.0 # DON'T DO THIS!!! THIS HIDES ERRORS
    res
end

# version that stores the simulated μ for each iteration
function ldpred_gibbs_acc(z, D, p, σ2, μ0; n_burnin = 100, n_iter = 500, verbose = false)

    n_variants = length(z)
    μ = μ0[:]
    Ω = zeros(Float64, n_variants)

    σ2_μ = σ2 / p
    σ_z = sqrt(1.0 + σ2_μ)
    η = σ2_μ / (1.0 + σ2_μ)  # shrinkage factor

    res = []
    
    for iter = 1:(n_burnin + n_iter)
        verbose && println("Iteration $iter/$(n_burnin + n_iter)")
        for j = 1:n_variants
            z̃_j = z[j] - μ'*D[:,j] + μ[j]
            p_j = 1.0 / (1.0 + (1.0-p)/p * σ_z * exp(-0.5 * z̃_j * z̃_j * η))
            μ[j] = (rand() <= p_j) ?  rand(Normal(η * z̃_j, sqrt(η))) : 0.0
        end

        if iter > n_burnin
            push!(res, copy(μ))
        end
    end


    res
end



#
# Empirical Bayes extension
#

# instead of using same prior for all μ_i, use prior as estimated from EB
#

function ldpred_gibbs_EB(z, D, π, σ_z, μ0, scale = 1.0; n_burnin = 100, n_iter = 500, verbose = false)

    K = length(π)
    n_variants = length(z)
    μ = μ0[:]
    Ω = zeros(Float64, n_variants)

    σ2_μ = max.(σ_z .^ 2 .- 1, 0.0)/scale
    σ_μ = sqrt.(σ2_μ)
    σ_z_scaled = sqrt.(1.0 .+ σ2_μ)
    η = σ2_μ ./ (1.0 .+ σ2_μ) # η[1,k] : shrinkage factor for component k
    
    for iter = 1:(n_burnin + n_iter)
        verbose && println("Iteration $iter/$(n_burnin + n_iter)")
        for j = 1:n_variants
            z̃_j = z[j] - μ'*D[:,j] + μ[j]
            p_j = estimate_component_prob([z̃_j], π, σ_z_scaled)
            
            k = sample(1:K, Weights(p_j[:]))
            μ[j] = rand(Normal(η[k] * z̃_j, sqrt(η[k])))
            #println("$j\t$(z[j])\t$(z̃_j)\t$(p_j)\t$(μ[j])")
            #if !isfinite(p_j[1])
            #    break
            #end
            if iter > n_burnin
                Ω[j] += sum(p_j .* η) * z̃_j 
            end
        end
    end

    res = Ω ./ n_iter
    res[(!isfinite).(res)] .= 0.0 # DON'T DO THIS!!! THIS HIDES ERRORS
    res
end


