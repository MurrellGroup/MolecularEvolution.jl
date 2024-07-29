
function univariate_modifier(f, modifier::UnivariateSampler; curr_value=nothing, kwargs...)
    return univariate_sampler(f, modifier, curr_value)
end

struct BranchlengthSampler <: UnivariateSampler
    #The first entry in acc_ratio holds the number of accepted proposals and the second entry holds the number of rejected proposals.
    acc_ratio
    log_bl_proposal
    log_bl_prior
    BranchlengthSampler(log_bl_proposal,log_bl_prior) = new([0,0],log_bl_proposal,log_bl_prior)
end 

"""  
    univariate_sampler(LL, modifier::BranchlengthPeturbation, curr_branchlength)

A MCMC algorithm that draws the next sample of a Markov Chain that approximates the Posterior distrubution over the branchlengths.
"""
function univariate_sampler(LL, modifier::BranchlengthSampler, curr_branchlength)
    return branchlength_metropolis(LL, modifier, curr_branchlength)
end

