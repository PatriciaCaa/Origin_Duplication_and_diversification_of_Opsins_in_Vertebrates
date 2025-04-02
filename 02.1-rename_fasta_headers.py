### RENAME FASTA HEADERS
## Project: Origin, duplication, and diversification of opnsin genes in vertebrates

# This script renames duplicate FASTA headers to ensure each sequence name is unique.

from collections import defaultdict

input_file = "AllOpsinhits80-melatonin.fasta.mafft-align-0.2.trimmed"
output_file = "AllOpsinhits80-melatonin.fasta.mafft-align-0.2.trimmed.renamed"

sequence_counts = defaultdict(int)

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            sequence_name = line.strip()
            sequence_counts[sequence_name] += 1
            if sequence_counts[sequence_name] > 1:
                sequence_name = f"{sequence_name}_{sequence_counts[sequence_name]}"
            outfile.write(sequence_name + "\n")
        else:
            outfile.write(line)

print(f"Renamed sequences have been saved to {output_file}")

