
#This should contain all code for going from sequence (or other?) data to partition vectors. This could be faster, but isn't currently a bottleneck.
#Need a nice way of specifying rules to cover all sensible ways of handling gap characters (where eg. gaps on the ends might have a different meaning)
#Need to have (chosing sensible conventions etc):
#0: Nuc data.
#1: Codon data.
#3: Amino acid data.
#4: Maybe date from leaf name functions?
#5: Country state from leaf name functions?
#6: Code to handle when the state comes in from a .csv with associated names, or something.

const nucstring = "ACGT"
const gappynucstring = "ACGT-"
const AAstring = "ACDEFGHIKLMNPQRSTVWY"
const gappyAAstring = "ACDEFGHIKLMNPQRSTVWY-"

#Dictionary messages default to Float64. You need to cast them to your required type in the use functions (eg. string2partition)
function one_of_k_dict_from_iter(st)
    dic = Dict{Char,Vector{Float64}}(Dict())
    empty = zeros(length(st))
    for (i,c) in enumerate(st)
        dic[c] = copy(empty)
        dic[c][i] = 1.0
    end
    return dic
end

AA_dict = one_of_k_dict_from_iter(AAstring)
gappy_AA_dict = one_of_k_dict_from_iter(gappyAAstring)
#Need to construct the IUPAC dictionaries customly.
nuc_dict = Dict{Char,Vector{Float64}}(Dict())
nuc_dict['A'] = [1.0, 0.0, 0.0, 0.0]
nuc_dict['C'] = [0.0, 1.0, 0.0, 0.0]
nuc_dict['G'] = [0.0, 0.0, 1.0, 0.0]
nuc_dict['T'] = [0.0, 0.0, 0.0, 1.0]
nuc_dict['B'] = [0.0, 1.0, 1.0, 1.0]
nuc_dict['M'] = [1.0, 1.0, 0.0, 0.0]
nuc_dict['Y'] = [0.0, 1.0, 0.0, 1.0]
nuc_dict['D'] = [1.0, 0.0, 1.0, 1.0]
nuc_dict['V'] = [1.0, 1.0, 1.0, 0.0]
nuc_dict['U'] = [0.0, 0.0, 0.0, 1.0]
nuc_dict['H'] = [1.0, 1.0, 0.0, 1.0]
nuc_dict['K'] = [0.0, 0.0, 1.0, 1.0]
nuc_dict['R'] = [1.0, 0.0, 1.0, 0.0]
nuc_dict['W'] = [1.0, 0.0, 0.0, 1.0]
nuc_dict['S'] = [0.0, 1.0, 1.0, 0.0]
nuc_dict['N'] = [1.0, 1.0, 1.0, 1.0]
nuc_dict['-'] = [1.0, 1.0, 1.0, 1.0]
gappy_nuc_dict = Dict{Char,Vector{Float64}}(Dict())
gappy_nuc_dict['A'] = [1.0, 0.0, 0.0, 0.0, 0.0]
gappy_nuc_dict['C'] = [0.0, 1.0, 0.0, 0.0, 0.0]
gappy_nuc_dict['G'] = [0.0, 0.0, 1.0, 0.0, 0.0]
gappy_nuc_dict['T'] = [0.0, 0.0, 0.0, 1.0, 0.0]
gappy_nuc_dict['B'] = [0.0, 1.0, 1.0, 1.0, 0.0]
gappy_nuc_dict['M'] = [1.0, 1.0, 0.0, 0.0, 0.0]
gappy_nuc_dict['Y'] = [0.0, 1.0, 0.0, 1.0, 0.0]
gappy_nuc_dict['D'] = [1.0, 0.0, 1.0, 1.0, 0.0]
gappy_nuc_dict['V'] = [1.0, 1.0, 1.0, 0.0, 0.0]
gappy_nuc_dict['U'] = [0.0, 0.0, 0.0, 1.0, 0.0]
gappy_nuc_dict['H'] = [1.0, 1.0, 0.0, 1.0, 0.0]
gappy_nuc_dict['K'] = [0.0, 0.0, 1.0, 1.0, 0.0]
gappy_nuc_dict['R'] = [1.0, 0.0, 1.0, 0.0, 0.0]
gappy_nuc_dict['W'] = [1.0, 0.0, 0.0, 1.0, 0.0]
gappy_nuc_dict['S'] = [0.0, 1.0, 1.0, 0.0, 0.0]
gappy_nuc_dict['N'] = [1.0, 1.0, 1.0, 1.0, 0.0]
gappy_nuc_dict['-'] = [0.0, 0.0, 0.0, 0.0, 1.0]

function extract(part::DiscretePartition,code::String)
    #NOTE: NEED TO HANDLE CASES WHERE TWO STATES ARE EQUAL. RETURN AMBIG CHAR? MAYBE CUSTOM BEHAVIOR FOR IUPAC HERE.
    code_arr = collect(code) #Need to compare this to just indexing into the string?
    return join([code_arr[argmax(part.state[:,i])] for i in 1:part.sites])
end
function extract(part::NucleotidePartition)
    extract(part,nucstring)
end
function extract(part::GappyNucleotidePartition)
    extract(part,gappynucstring)
end
function extract(part::AminoAcidPartition)
    extract(part,AAstring)
end
function extract(part::GappyAminoAcidPartition)
    extract(part,gappyAAstring)
end



#Creates new memory. This is currently required for the populate_tree function to do things automatically.
function obs2partition!(dest::DiscretePartition,seq::String,dic::Dict)
    #if length(seq) != dest.sites
    #    @warn "obs2partition!() is changing the number of sites in a partition."
    #end
    dest.state = zeros(eltype(dest.state),dest.states,length(seq))
    dest.sites = length(seq)
    dest.scaling = zeros(length(seq))
    unmatched = ones(eltype(dest.state),dest.states)
    for (i,c) in enumerate(uppercase(seq))
        dest.state[:,i] = get(dic,c,unmatched)
    end
end
function obs2partition!(dest::NucleotidePartition,seq::String)
    obs2partition!(dest,seq,nuc_dict)
end
function obs2partition!(dest::GappyNucleotidePartition,seq::String)
    obs2partition!(dest,seq,gappy_nuc_dict)
end
function obs2partition!(dest::AminoAcidPartition,seq::String)
    obs2partition!(dest,seq,AA_dict)
end
function obs2partition!(dest::GappyAminoAcidPartition,seq::String)
    obs2partition!(dest,seq,gappy_AA_dict)
end

