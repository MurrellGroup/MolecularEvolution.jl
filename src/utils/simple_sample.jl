
function univariate_modifier(f, modifier::UnivariateSampler; curr_value=nothing, kwargs...)
    return univariate_sampler(f, modifier, curr_value)
end

"""
    BranchlengthSampler

    A type that allows you to specify a additive proposal function in the log domain and a prior distrubution over the log of the branchlengths. It also stores `acc_ratio` which is a tuple of `(ratio, total, #acceptances)`, where `ratio::Float64` is the acceptance ratio, `total::Int64` is the total number of proposals, and `#acceptances::Int64` is the number of acceptances.
"""
mutable struct BranchlengthSampler <: UnivariateSampler
    acc_ratio::Tuple{Float64, Int64, Int64} #(ratio, total, #acceptances)
    log_bl_proposal::Distribution
    log_bl_prior::Distribution
    BranchlengthSampler(log_bl_proposal,log_bl_prior) = new((0.0, 0, 0),log_bl_proposal,log_bl_prior)
end 

"""  
    univariate_sampler(LL, modifier::BranchlengthPeturbation, curr_branchlength)

    A MCMC algorithm that draws the next sample of a Markov Chain that approximates the Posterior distrubution over the branchlengths.
"""
function univariate_sampler(LL, modifier::BranchlengthSampler, curr_branchlength)
    return metropolis_step(LL, modifier, curr_branchlength)
end

"""
    metropolis_step(LL::Function, modifier, curr_value)

Does a standard metropolis step in an MCMC, i.e. proposes a candidate symmetrically and returns the next
state in the chain, decided by the candidate being rejected or not.
# Interface
You need a `MySampler <: Any` to implement
- `proposal(modifier::MySampler, curr_value)`
- `log_prior(modifier::MySampler, x)`
- `apply_decision(modifier::MySampler, accept::Bool)`

`LL` is by default called on `curr_value` and the returned value of `proposal`.
Although, it is possible to transform the current value before proposing a new value, and
then take the inverse transform to match the argument `LL` expects.
# Extended interface
To make proposals in a transformed space, you overload
- `tr(modifier::MySampler, x)`
- `invtr(modifier::MySampler, x)`
which are identity transformations by default.
"""
function metropolis_step(LL::Function, modifier, curr_value)
    tr_curr_value = tr(modifier, curr_value)
    tr_prop = proposal(modifier, tr_curr_value)
    accept_proposal =
        rand() <= exp(
            LL(invtr(modifier, tr_prop)) + log_prior(modifier, tr_prop) -
            LL(invtr(modifier, tr_curr_value)) - log_prior(modifier, tr_curr_value),
        )
    apply_decision(modifier, accept_proposal)
    return invtr(modifier, ifelse(accept_proposal, tr_prop, tr_curr_value))
end

# The prior distribution for the variable log(branchlength). A small perturbation of +1e-12 is added to enhance numerical stability near zero.
proposal(modifier::BranchlengthSampler, curr_value) = exp(log(curr_value) + rand(modifier.log_bl_proposal))
log_prior(modifier::BranchlengthSampler, x) = logpdf(modifier.log_bl_prior, log(x + 1e-12))
#Default definition. Overload it for your own modifier type
function apply_decision(modifier, accept::Bool)
    ratio, total, acc = modifier.acc_ratio
    total += 1
    if accept
        acc += 1
    end
    ratio = acc / total
    modifier.acc_ratio = (ratio, total, acc)
end

tr(modifier, x) = x
invtr(modifier, x) = x
