import .FASTX

export read_fasta
"""
    read_fasta(filepath::String)

Reads in a fasta file and returns a tuple of (seqnames, seqs).
"""
function read_fasta(filepath::String)
    reader = FASTX.FASTA.Reader(open(filepath, "r"))
    fasta_in = [record for record in reader]
    close(reader)
    return [String(FASTX.FASTA.identifier(rec)) for rec in fasta_in],
    [uppercase(String(FASTX.FASTA.sequence(rec))) for rec in fasta_in]
end

export write_fasta
"""
    write_fasta(filepath::String, sequences::Vector{String}; seq_names = nothing)

Writes a fasta file from a vector of sequences, with optional seq_names.
"""
function write_fasta(filepath::String, sequences::Vector{String}; seq_names = nothing)
    if seq_names === nothing
        seq_names = ["S$(i)" for i = 1:length(sequences)]
    end
    writer = FASTX.FASTA.Writer(open(filepath, "w"))
    for i = 1:length(seq_names)
        rec = FASTX.FASTA.Record(seq_names[i], sequences[i])
        write(writer, rec)
    end
    close(writer)
end
