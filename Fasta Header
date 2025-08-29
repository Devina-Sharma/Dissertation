from Bio import SeqIO 

# Input and output file names
input_fasta = "rbcl_original.fasta"
output_fasta = "rbcl_clean.fasta"

# Open the output file to write the cleaned sequences
with open(output_fasta, "w") as out_handle:
       for record in SeqIO.parse(input_fasta, "fasta"):
                header = record.description

        try:
            # Try to get the part of the header after 'tax='
            tax_fields = header.split("tax=")[1]

            # Look through each tax field (separated by commas)
            for field in tax_fields.split(","):
                # If the field starts with 's:', it's the species name
                if field.startswith("s:"):
                    # Save the species name (after 's:')
                    species = field.split(":")[1]
                    break
            else:
                # If 's:' not found, use 'unknown_species'
                species = "unknown_species"

        except IndexError:
            # If there's no 'tax=' in the header, use 'unknown_species'
            species = "unknown_species"

        # Set the new name to the species (in lowercase)
        record.id = species.lower()

        # Remove the full description 
        record.description = ""

        # Write the cleaned species to the output file
        SeqIO.write(record, out_handle, "fasta")
