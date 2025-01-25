nuc2num = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4);

#Returns the "codon_pos,from,to" position in a nucleotide table
function codon_diff(c1, c2)
    diffs = [c1[i] != c2[i] for i = 1:3]
    if sum(diffs) == 1
        ind = [1, 2, 3][diffs][1]
        return (ind, nuc2num[c1[ind]], nuc2num[c2[ind]])
    else
        return (-1, -1, -1) #Ugh, gross.
    end
end

#A struct to hold the contents/lookup tables for a particular genetic code.
#Need to re-think what we should store here for various use-cases, including:
#Construction of codon Q matrices
#Construction of single rows of codon Q matrices
#Convert from codon to AA
#Convert from codon to nuc string
#F3x4 calculations
#complete exhaustive list
struct GeneticCode
    genetic_code::Dict{String,Char}
    codons::Array{String,1}
    sense_codons::Array{String,1}
    stop_codons::Array{String,1}
    amino_acid_lookup::Array{Char,1}
    sense2all::Array{Int64,1}
    string2sense::Dict{String,Int}
    syn_positions::Array{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64,Int64}},1}
    nonsyn_positions::Array{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64,Int64}},1}
    amino_acids::Array{Char,1}
    codon2AA_pos::Array{Int64,1}
    codon_to_nuc_map::Array{Tuple{Int64,Int64,Int64},2}
    num_sense::Int64

    function GeneticCode(genetic_code::Dict{String,Char})
        codons = sort(collect(keys(genetic_code)))
        sense_codons = [i for i in codons if genetic_code[i] != '*']
        stop_codons = [i for i in codons if genetic_code[i] == '*']
        amino_acid_lookup = [genetic_code[i] for i in sense_codons]
        sense2all = zeros(Int, length(sense_codons))
        count = 1
        for i = 1:length(codons)
            sense2all[count] = i
            if genetic_code[codons[i]] != '*'
                count += 1
            end
        end
        string2sense = Dict{String,Int}()
        for cod = 1:length(sense_codons)
            string2sense[sense_codons[cod]] = cod
        end

        #These encode the position in the codon model, the 1,2,3 codon position difference, and the corresponding nucleotide change.
        #They could be structured by nuc change to enable partial matrix updates when only some nuc parameters change.
        syn_positions = Array{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64,Int64}},1}([])
        nonsyn_positions = Array{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64,Int64}},1}([])
        #A dense version of the above, for using when you only need to construct one row at a time
        codon_to_nuc_map = hcat(
            [
                [(-1, -1, -1) for i = 1:length(sense_codons)] for
                i = 1:length(sense_codons)
            ]...,
        )

        for i = 1:length(sense_codons)
            for j = 1:length(sense_codons)
                c1, c2 = codons[sense2all[i]], codons[sense2all[j]]
                diff = codon_diff(c1, c2)
                if diff != (-1, -1, -1)

                    codon_to_nuc_map[i, j] = diff

                    if amino_acid_lookup[i] == amino_acid_lookup[j]
                        push!(syn_positions, ((i, j), diff))
                    else
                        push!(nonsyn_positions, ((i, j), diff))
                    end
                end
            end
        end

        amino_acids = sort(union([genetic_code[cod] for cod in sense_codons]))

        codon2AA_pos =
            [findfirst(genetic_code[cod] .== amino_acids) for cod in sense_codons]

        new(
            genetic_code,
            codons,
            sense_codons,
            stop_codons,
            amino_acid_lookup,
            sense2all,
            string2sense,
            syn_positions,
            nonsyn_positions,
            amino_acids,
            codon2AA_pos,
            codon_to_nuc_map,
            length(sense_codons),
        )
    end
end

universal_genetic_code = Dict(
    "ATA" => 'I',
    "ATC" => 'I',
    "ATT" => 'I',
    "ATG" => 'M',
    "ACA" => 'T',
    "ACC" => 'T',
    "ACG" => 'T',
    "ACT" => 'T',
    "AAC" => 'N',
    "AAT" => 'N',
    "AAA" => 'K',
    "AAG" => 'K',
    "AGC" => 'S',
    "AGT" => 'S',
    "AGA" => 'R',
    "AGG" => 'R',
    "CTA" => 'L',
    "CTC" => 'L',
    "CTG" => 'L',
    "CTT" => 'L',
    "CCA" => 'P',
    "CCC" => 'P',
    "CCG" => 'P',
    "CCT" => 'P',
    "CAC" => 'H',
    "CAT" => 'H',
    "CAA" => 'Q',
    "CAG" => 'Q',
    "CGA" => 'R',
    "CGC" => 'R',
    "CGG" => 'R',
    "CGT" => 'R',
    "GTA" => 'V',
    "GTC" => 'V',
    "GTG" => 'V',
    "GTT" => 'V',
    "GCA" => 'A',
    "GCC" => 'A',
    "GCG" => 'A',
    "GCT" => 'A',
    "GAC" => 'D',
    "GAT" => 'D',
    "GAA" => 'E',
    "GAG" => 'E',
    "GGA" => 'G',
    "GGC" => 'G',
    "GGG" => 'G',
    "GGT" => 'G',
    "TCA" => 'S',
    "TCC" => 'S',
    "TCG" => 'S',
    "TCT" => 'S',
    "TTC" => 'F',
    "TTT" => 'F',
    "TTA" => 'L',
    "TTG" => 'L',
    "TAC" => 'Y',
    "TAT" => 'Y',
    "TAA" => '*',
    "TAG" => '*',
    "TGC" => 'C',
    "TGT" => 'C',
    "TGA" => '*',
    "TGG" => 'W',
)

const universal_code = GeneticCode(universal_genetic_code);

function count_F3x4(seqs::Array{String})
    F3x4 = zeros(3, 4)
    for k = 1:3
        pos1inds = [3 * (i - 1) + k for i = 1:Int(length(seqs[1]) / 3)]
        d = proportionmap(collect(join([s[pos1inds] for s in seqs])))
        F3x4[k, 1] = get(d, 'A', 0.0)
        F3x4[k, 2] = get(d, 'C', 0.0)
        F3x4[k, 3] = get(d, 'G', 0.0)
        F3x4[k, 4] = get(d, 'T', 0.0)
    end
    return F3x4
end

function F3x4_eq_freqs(F3x4; genetic_code = universal_code)
    eq = [
        F3x4[1, nuc2num[c[1]]] * F3x4[2, nuc2num[c[2]]] * F3x4[3, nuc2num[c[3]]] for
        c in genetic_code.sense_codons
    ]
    return eq ./ sum(eq)
end

function MG94_F3x4(alpha, beta, nuc_matrix, F3x4; genetic_code = universal_code)
    codon_matrix =
        zeros(length(genetic_code.sense_codons), length(genetic_code.sense_codons))
    for p in genetic_code.syn_positions
        codon_matrix[p[1][1], p[1][2]] =
            alpha * nuc_matrix[p[2][2], p[2][3]] * F3x4[p[2][1], p[2][3]]
    end
    for p in genetic_code.nonsyn_positions
        codon_matrix[p[1][1], p[1][2]] =
            beta * nuc_matrix[p[2][2], p[2][3]] * F3x4[p[2][1], p[2][3]]
    end
    for i = 1:length(genetic_code.sense_codons)
        codon_matrix[i, i] = -sum(codon_matrix[i, :])
    end
    return codon_matrix
end

function count_F61(seqs::Array{String}; code = universal_code)
    F61 = zeros(length(code.sense_codons))
    dic = proportionmap(
        vcat(
            [[seq[3*(i-1)+1:3*(i-1)+3] for i = 1:Int(length(seq) / 3)] for seq in seqs]...,
        ),
    )
    for i = 1:length(F61)
        F61[i] = get(dic, code.sense_codons[i], 0.0)
    end
    return sum2one(F61)
end

function MG94_F61(alpha, beta, nuc_matrix, F61; genetic_code = universal_code)
    codon_matrix =
        zeros(length(genetic_code.sense_codons), length(genetic_code.sense_codons))
    for p in genetic_code.syn_positions
        codon_matrix[p[1][1], p[1][2]] = alpha * nuc_matrix[p[2][2], p[2][3]] * F61[p[1][2]]
    end
    for p in genetic_code.nonsyn_positions
        codon_matrix[p[1][1], p[1][2]] = beta * nuc_matrix[p[2][2], p[2][3]] * F61[p[1][2]]
    end
    for i = 1:length(genetic_code.sense_codons)
        codon_matrix[i, i] = -sum(codon_matrix[i, :])
    end
    return codon_matrix
end

function HB98_F61(alpha, nuc_matrix, F61; genetic_code = universal_code)
    #See https://www.ncbi.nlm.nih.gov/pubmed/9656490 but beware of the typo in equation 9.
    codon_matrix =
        zeros(length(genetic_code.sense_codons), length(genetic_code.sense_codons))

    for p in genetic_code.syn_positions
        #a substitution from codon p[1][1] to codon p[1][2]

        PIaPab = F61[p[1][1]] * nuc_matrix[p[2][2], p[2][3]]
        PIbPba = F61[p[1][2]] * nuc_matrix[p[2][3], p[2][2]]
        if PIaPab != PIbPba
            f_ab = log(PIbPba / PIaPab) / (1 - (PIaPab / PIbPba))
        else
            f_ab = 1.0
        end
        codon_matrix[p[1][1], p[1][2]] = alpha * nuc_matrix[p[2][2], p[2][3]] * f_ab
    end
    for p in genetic_code.nonsyn_positions
        #duplicated code, because this model doesn't really use the syn/non-syn distinction
        PIaPab = F61[p[1][1]] * nuc_matrix[p[2][2], p[2][3]]
        PIbPba = F61[p[1][2]] * nuc_matrix[p[2][3], p[2][2]]
        if PIaPab != PIbPba
            f_ab = log(PIbPba / PIaPab) / (1 - (PIaPab / PIbPba))
        else
            f_ab = 1.0
        end
        codon_matrix[p[1][1], p[1][2]] = alpha * nuc_matrix[p[2][2], p[2][3]] * f_ab
    end
    for i = 1:length(genetic_code.sense_codons)
        codon_matrix[i, i] = -sum(codon_matrix[i, :])
    end
    return codon_matrix
end

function HB98_AAfit(alpha, nuc_matrix, AA_fitness; genetic_code = universal_code)
    #Halpern and Bruno, but with the fitnesses directly parameterized.
    #Closer to https://academic.oup.com/mbe/article/32/4/1097/1077799
    #2Ns_ij = 2Nf_j - 2Nf_i
    #Sub rate = mutation_ij * 2Ns_ij/(1-e^-2Ns_ij)
    codon_matrix =
        zeros(length(genetic_code.sense_codons), length(genetic_code.sense_codons))
    for p in genetic_code.syn_positions
        #In this AA_fitness model, fitnesses are always equal for syn changes.
        codon_matrix[p[1][1], p[1][2]] = alpha * nuc_matrix[p[2][2], p[2][3]]
    end
    for p in genetic_code.nonsyn_positions
        #duplicated code, because this model doesn't really use the syn/non-syn distinction
        diff = (
            AA_fitness[genetic_code.codon2AA_pos[p[1][2]]] -
            AA_fitness[genetic_code.codon2AA_pos[p[1][1]]]
        )#/(1+AA_fitness[genetic_code.codon2AA_pos[p[1][1]]]) #Need to think about this term in the context of the approximation
        if abs(diff) < 0.001 #To catch the hole discontinuity a 2nd order approximation that is REALLY close.
            f_ab = 1 / 12 * (12 + diff * (6 + diff)) #1+0.5*diff
        else
            f_ab = diff / (1 - exp(-diff))
        end
        codon_matrix[p[1][1], p[1][2]] = alpha * nuc_matrix[p[2][2], p[2][3]] * f_ab
    end
    for i = 1:length(genetic_code.sense_codons)
        codon_matrix[i, i] = -sum(codon_matrix[i, :])
    end
    return codon_matrix
end

function codonHKY85(TrTv)
    return [
        0.0 1.0 TrTv 1.0
        1.0 0.0 1.0 TrTv
        TrTv 1.0 0.0 1.0
        1.0 TrTv 1.0 0.0
    ]
end

mutable struct CodonPartition <: DiscretePartition
    state::Array{Float64,2}
    states::Int
    sites::Int
    scaling::Array{Float64,1}
    function CodonPartition(sites; code = universal_code)
        new(
            zeros(length(code.sense_codons), sites),
            length(code.sense_codons),
            sites,
            zeros(sites),
        )
    end
    function CodonPartition(state, states, sites, scaling; code = universal_code)
        @assert size(state) == (states, sites) && states == length(code.sense_codons)
        new(state, states, sites, scaling)
    end
end

#Make this handle IUPAC ambigs sensible. Any codon compatible with the ambig should get a 1.0
function obs2partition!(dest::CodonPartition, seq::String; code = universal_code)
    problem_codons = String[]
    if mod(length(seq), 3) != 0
        error("Codon sequences must be divisible by 3")
    end
    if length(seq) / 3 != dest.sites
        error("Sequence length does not match partition")
    end

    @views for j in axes(dest.state, 2)
        c = seq[3*(j-1)+1:3*(j-1)+3]
        cod_ind = get(code.string2sense, c, -1)
        if cod_ind == -1
            fill!(dest.state[:, j], 1.0)
            push!(problem_codons, c)
        else
            fill!(dest.state[:, j], 0.0)
            dest.state[cod_ind, j] = 1.0
        end
    end
    fill!(dest.scaling, 0.0)
    return countmap(problem_codons)
end


function partition2obs(part::CodonPartition; code = universal_code)
    return join([code.sense_codons[argmax(part.state[:, i])] for i = 1:part.sites])
end
