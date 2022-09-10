# Optimization

There are two distinct kinds of optimization: "global" model parameters, and then tree branchlengths and topology. These are kept distinct because we can use algorithmic tricks to dramatically improve the performance of the latter.

## Optimizing model parameters

We first need to construct an objective function. A very common use case involves parameterizing a rate matrix (along with all the constraints this entails) from a flat parameter vector. `reversibleQ` can be convenient here, which takes a vector of parameters and equilibrium frequencies and returns a reversible rate matrix.
The parameters are the upper triangle (excluding the diagonal) of the rate matrix:

```@example 1
using MolecularEvolution #hide
reversibleQ(1:6,ones(4))
```

...and the equilibrium frequencies are multiplied column-wise:

```@example 1
reversibleQ(ones(6),[0.1,0.2,0.3,0.4])
```

Another convenient trick is to be able to parameterize a vector of positive frequencies that sum to 1, using N-1 unconstrained parameters. `unc2probvec` can help:

```@example 1
unc2probvec(zeros(3))
```

[`ParameterHandling.jl`](https://github.com/invenia/ParameterHandling.jl) provides a convenient framework for managing collections of parameters in a way that plays with much of the Julia optimization ecosystem, and we recommend its use. Here we'll use `ParameterHandling` and [`NLopt`](https://github.com/JuliaOpt/NLopt.jl).

First, we'll load in some example nucleotide data:

```julia
using MolecularEvolution, FASTX, ParameterHandling, NLopt

#Read in seqs and tree, and populate the three  NucleotidePartitions
seqnames, seqs = read_fasta("Data/MusNuc_IGHV.fasta")
tree = read_newick_tree("Data/MusNuc_IGHV.tre")
initial_partition = NucleotidePartition(length(seqs[1]))
populate_tree!(tree,initial_partition,seqnames,seqs)
```

Then we set up the model parameters, and the objective function:

```julia
#Named tuple of parameters, with initial values and constraints (from ParameterHandling.jl)
initial_params = (
        rates=positive(ones(6)), #rates must be non-negative
        pi=zeros(3) #will be transformed into 4 eq freqs
)
flat_initial_params, unflatten = value_flatten(initial_params) #See ParameterHandling.jl docs

#Set up a function that builds a model from these parameters
function build_model_vec(params)
    pi = unc2probvec(params.pi)
    return DiagonalizedCTMC(reversibleQ(params.rates,pi))
end

#Set up the function to be *minimized*
function objective(params::NamedTuple; tree = tree)
    #In this example, we are optimizing the nuc equilibrium freqs
    #We'll also assume that the starting frequencies (at the root of the tree) are the eq freqs
    tree.parent_message[1].state .= unc2probvec(params.pi)
    return -log_likelihood!(tree,build_model_vec(params)) #Note, negative of LL, because minimization
end

num_params = length(flat_initial_params)
```

Then we'll set up an optimizer from `NLOpt`. See [this discussion](https://discourse.julialang.org/t/optim-what-optimiser-is-best-if-your-gradient-computation-is-slow/5487/12) and [this exploration](https://github.com/SciML/DiffEqParamEstim.jl/blob/master/test/lorenz_true_test.jl) of optimizers.

```julia
#See: 
opt = Opt(:LN_BOBYQA, num_params)
#Note: NLopt requires a function that returns a gradient, even for gradient free methods, hence (x,y)->...
min_objective!(opt, (x,y) -> (objective ∘ unflatten)(x)) #See ParameterHandling.jl docs for objective ∘ unflatten explanation
#Some bounds (which will be in the transformed domain) to prevent searching numerically silly bits of parameter space:
lower_bounds!(opt, [-10.0 for i in 1:num_params])
upper_bounds!(opt, [10.0 for i in 1:num_params])
xtol_rel!(opt, 1e-12)
_,mini,_ = NLopt.optimize(opt, flat_initial_params)
final_params = unflatten(mini)

optimzed_model = build_model_vec(final_params)
println("Opt LL:",log_likelihood!(tree,optimzed_model))
```

We can view the optimized parameter values:

```
@show final_params.rates
```

Or the entire optimized rate matrix:

```julia
matrix_for_display(optimzed_model.Q,['A','C','G','T'])
```
```
Opt LL:-3783.226756522292
5×5 Matrix{Any}:
 ""     'A'        'C'        'G'        'T'
 'A'  -1.02672    0.246386   0.494024   0.286309
 'C'   0.314289  -0.971998   0.23034    0.427368
 'G'   0.587774   0.214842  -0.950007   0.147391
 'T'   0.300663   0.35183    0.130093  -0.782586
```



