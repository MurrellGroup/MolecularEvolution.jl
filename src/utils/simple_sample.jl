
function univariate_modifier(f, modifier::UnivariateSampler; curr_value=nothing, kwargs...)
    return univariate_sampler(f, modifier, curr_value)
end

struct BranchlengthPerturbation <: UnivariateSampler
    sigma
    #The first entry in acc_ratio holds the number of accepted proposals and the second entry holds the number of rejected proposals.
    acc_ratio
    BranchlengthPerturbation(sigma) = new(sigma, [0,0])
end 

"""  
    univariate_sampler(LL, modifier::BranchlengthPeturbation, curr_branchlength)

A MCMC algorithm that draws the next sample of a Markov Chain that approximates the Posterior distrubution over the branchlengths.
"""
function univariate_sampler(LL, modifier::BranchlengthPerturbation, curr_branchlength)
    return branchlength_metropolis(LL, modifier, curr_branchlength)
end

