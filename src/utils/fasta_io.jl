import .FASTX

export read_fasta
function read_fasta(filepath::String)
    reader = FASTX.FASTA.Reader(open(filepath, "r"))
    fasta_in = [record for record in reader]
    close(reader)
    return [String(FASTX.FASTA.identifier(rec)) for rec in fasta_in],
    [String(FASTX.FASTA.sequence(rec)) for rec in fasta_in]
end

export write_fasta
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
