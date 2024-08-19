
function univariate_modifier(f, modifier::UnivariateSampler; curr_value=nothing, kwargs...)
    return univariate_sampler(f, modifier, curr_value)
end
"""  
    univariate_sampler(LL, modifier::BranchlengthPeturbation, curr_branchlength)

A MCMC algorithm that draws the next sample of a Markov Chain that approximates the Posterior distrubution over the branchlengths.
"""
function univariate_sampler(LL, modifier::BranchlengthSampler, curr_branchlength)
    return branchlength_metropolis(LL, modifier, curr_branchlength)
end

struct BranchlengthSampler <: UnivariateSampler
    #The first entry in acc_ratio holds the number of accepted proposals and the second entry holds the number of rejected proposals.
    acc_ratio
    log_bl_proposal
    log_bl_prior
    BranchlengthSampler(log_bl_proposal,log_bl_prior) = new([0,0],log_bl_proposal,log_bl_prior)
end 

function branchlength_metropolis(LL, modifier, curr_value)
    # The prior distribution for the variable log(branchlength). A small perturbation of +1e-12 is added to enhance numerical stability near zero.
    log_prior(x) = logpdf(modifier.log_bl_prior,log(x + 1e-12))
    # Adding additive normal symmetrical noise in the log(branchlength) domain to ensure the proposal function is symmetric.
    noise = rand(modifier.log_bl_proposal)
    proposal = exp(log(curr_value)+noise)
    # The standard Metropolis acceptance criterion.
    if rand() <= exp(LL(proposal)+log_prior(proposal)-LL(curr_value)-log_prior(curr_value))
        modifier.acc_ratio[1] = modifier.acc_ratio[1] + 1
        return proposal
    else
        modifier.acc_ratio[2] = modifier.acc_ratio[2] + 1 
        return curr_value
    end
end




