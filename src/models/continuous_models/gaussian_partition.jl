#-----------------Gaussian prop---------------
#Partition behavior
mutable struct GaussianPartition <: ContinuousPartition
    mean::Float64
    var::Float64
    norm_const::Float64
    function GaussianPartition(mean,var,norm_const)
        new(mean,var,norm_const)
    end
    function GaussianPartition(mean,var)
        new(mean,var,0.0)
    end
    function GaussianPartition()
        new(0.0,1.0,0.0)
    end
end

#From the first section of http://www.tina-vision.net/docs/memos/2003-003.pdf
function merge_two_gaussians(g1::GaussianPartition, g2::GaussianPartition)
    #Handling some edge cases. These aren't mathematically sensible. A gaussian with "Inf" variance will behave like a 1,1,1,1 vector in discrete felsenstein.
    #To-do: update some of these so that the norm constant is properly handled, even if the variance is Inf (so it isn't exactly well-defined anyway)
    if g1.var == 0 && g2.var == 0 && g1.var != g2.var
        error("both gaussians have 0 variance but different means")
    elseif g1.var == 0
        return deepcopy(g1)
    elseif g2.var == 0
        return deepcopy(g2)
    end
    if g1.var == Inf && g2.var == Inf
        return GaussianPartition((g1.mean+g2.mean)/2,Inf,0.0)
    elseif g1.var == Inf
        return deepcopy(g2)
    elseif g2.var == Inf
        return deepcopy(g1)
    end
    res_gaussian = GaussianPartition()
    res_gaussian.var = 1/(1/g1.var + 1/g2.var)
    res_gaussian.mean = res_gaussian.var * (g1.mean/g1.var + g2.mean/g2.var)
    # log of scaling constant
    res_gaussian.norm_const = -0.5 * ( log(2 * pi * (g1.var * g2.var / res_gaussian.var)) + (g1.mean^2 / g1.var) + (g2.mean^2 / g2.var) - (res_gaussian.mean^2 / res_gaussian.var) )
    res_gaussian.norm_const += (g1.norm_const + g2.norm_const)
    return res_gaussian
end

function identity!(dest::GaussianPartition)
    dest.var = Inf
    dest.norm_const = 0.0
end

function combine!(dest::GaussianPartition,src::GaussianPartition)
    new_g = merge_two_gaussians(dest,src)
    dest.mean = new_g.mean
    dest.var = new_g.var
    dest.norm_const = new_g.norm_const
end

function site_LLs(part::GaussianPartition)
    return [part.norm_const]
end

function gaussian_pdf(g::GaussianPartition, x::Float64)
    if g.var == 0
        return Float64(x == g.mean) #Hokey...
    end
    return pdf(Normal(g.mean, sqrt(g.var)), x)
end

#And sampling
function sample_partition!(partition::GaussianPartition)
    partition.mean = randn()*sqrt(partition.var) + partition.mean
    partition.var = 0.0
    partition.norm_const = 0.0
end

#And max
function max_partition!(partition::GaussianPartition)
    partition.var = 0.0
    partition.norm_const = 0.0
end
 
function obs2partition!(partition::GaussianPartition, obs::Float64)
    partition.mean = obs
    partition.var = 0.0
    partition.norm_const = 0.0
end

function extract(partition::GaussianPartition)
    return partition.mean
end