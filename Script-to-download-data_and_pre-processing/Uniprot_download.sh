### Script to download genomes from the UniProt database

#!/bin/bash

# Input file containing accession numbers
accession_file="Uniprot_accession_numbers.txt"

# Output directory for downloaded sequences
genome_sequences="uniprot_sequences"
mkdir -p "$genome_sequences"

# Uniprot base url
url="https://rest.uniprot.org/uniprotkb/"

# Loop through each accession number in the file
while read -r accession; do
    if [[ -n "$accession" ]]; then
        # Download the sequence in FASTA format
        echo "Downloading sequence for accession: $accession"
        curl -s "${url}${accession}.fasta" -o "${genome_sequences}/${accession}.fasta"
    fi
done < "$accession_file"

echo "All downloads complete. Sequences saved in: $genome_sequences"

