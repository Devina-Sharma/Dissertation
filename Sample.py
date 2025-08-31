import os
import random
from Bio import SeqIO  

# input FASTA file and output folder for simulated samples
input_fasta = "rbcl_trim.fasta"
output_folder = "rbcl_simulated_samples"

num_samples = 40               # Total number of simulated samples to generate
min_species = 1                # Minimum number of species per sample
max_species = 5                # Maximum number of species per sample
total_reads_per_sample = 80000  # Total number of reads per sample (paired-end reads)

# Create output directory
os.makedirs(output_folder, exist_ok=True)

# Load all sequences from the input FASTA file into a dictionary 
records = list(SeqIO.parse(input_fasta, "fasta"))
sequence_dict = {i: record for i, record in enumerate(records)}

print(f"Total sequences loaded: {len(sequence_dict)}")

# Generate simulated samples
for sample_num in range(1, num_samples + 1):
    
    # Randomly choose how many species to include in this sample
    num_species = random.randint(min_species, max_species)
    keepgoing = True
    while keepgoing:
        selected_ids = random.sample(list(sequence_dict.keys()), num_species)
        selected_seqs = [sequence_dict[i] for i in selected_ids]
        
        # Extract sequence IDs to check for duplicates
        seq_list = [seq.id for seq in selected_seqs]
        print(seq_list)  # Optional: print selected sequence IDs
        print(len(seq_list))
        print(len(set(seq_list)))
        
        # Continue if all selected sequence IDs are unique
        if len(seq_list) == len(set(seq_list)):
            keepgoing = False
        else:
            keepgoing = True

    # Write the selected sequences to a sample-specific FASTA file
    fasta_filename = os.path.join(output_folder, f"sample{sample_num}.fasta")
    with open(fasta_filename, "w") as fasta_file:
        for record in selected_seqs:
            fasta_file.write(f">{record.id}\n{record.seq}\n")

    # Assign random weight to each species 
    counts = [random.randint(1, 100) for _ in selected_seqs]
    total_counts = sum(counts)

    # Make sure each count is even (for paired-end reads)
    read_counts = [int((total_reads_per_sample * c / total_counts) // 2) * 2 for c in counts]

    # Adjust the last count to ensure total is exactly 80,000
    read_counts[-1] += total_reads_per_sample - sum(read_counts)

    count_filename = os.path.join(output_folder, f"sample{sample_num}_counts.txt")
    with open(count_filename, "w") as count_file:
        for record, count in zip(selected_seqs, read_counts):
            count_file.write(f"{record.id}\t{count}\n")


print("FASTA + count files created for all 40 samples with counts summing to 80,000.")
