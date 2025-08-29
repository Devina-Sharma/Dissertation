import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Input/output files
input_fasta = "rbcl_database.fasta"
filtered_fasta = "rbcl_filtered.fasta"
trimmed_fasta = "rbcl_trim.fasta"

# Primers (5' to 3')
forward_primer = "ATGTCACCACAAACAGAGACTAAAGC"
reverse_primer = "AGGGGACGACCATACTTGTTCA"

# IUPAC dictionary
iupac_dict = 
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

def iupac_to_regex(primer):
    return ''.join(iupac_dict.get(base, base) for base in primer.upper())

# Compile regex patterns
fwd_pattern = re.compile(iupac_to_regex(forward_primer))
rev_pattern = re.compile(iupac_to_regex(str(Seq(reverse_primer).reverse_complement())))

# Step 1: Filter sequences containing both primers
sequences_with_primers = []
for record in SeqIO.parse(input_fasta, "fasta"):
    seq_str = str(record.seq)

    if fwd_pattern.search(seq_str) and rev_pattern.search(seq_str):
            new_record = SeqRecord(
            record.seq,
            id=record.id,
            name=record.id,
            description=""
        )
        sequences_with_primers.append(new_record)
    else:
        print(f"Skipping {record.id}: primers not found.")

SeqIO.write(sequences_with_primers, filtered_fasta, "fasta")
print(f"\nSaved {len(sequences_with_primers)} sequences with primers to '{filtered_fasta}'.")

# Step 2: Trim between primers (include primers, remove extra)
trimmed_records = []
for record in sequences_with_primers:
    seq_str = str(record.seq)

    fwd_match = fwd_pattern.search(seq_str)
    rev_match = rev_pattern.search(seq_str)

    if fwd_match and rev_match:
        start = fwd_match.start()
        end = rev_match.end()

        if start < end:
            trimmed_seq = seq_str[start:end]
            trimmed_record = SeqRecord(
                Seq(trimmed_seq),
                id=record.id,
                name=record.id,
                description=""
            )
            trimmed_records.append(trimmed_record)
        else:
            print(f"Skipping {record.id}: forward primer found after reverse primer.")
    else:
        print(f"Skipping {record.id}: primer match failed during trimming.")

# Save trimmed sequences
if trimmed_records:
    SeqIO.write(trimmed_records, trimmed_fasta, "fasta")
    print(f"\nTrimmed {len(trimmed_records)} sequences written to '{trimmed_fasta}'.")
else:
    print("\nNo trimmed sequences created.")
