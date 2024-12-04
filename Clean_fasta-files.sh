### SCRIPT TO CLEAN FASTA FILES

# This script cleans up FASTA files (.faa) by removing unwanted characters from sequence lines.
# - Removes periods (.) that may represent gaps in alignments.
# - Removes asterisks (*) that may represent stop codons in protein sequences.
# This ensures the sequences are clean and compatible with downstream analyses.


for i in *faa 
do 
    sed -i '/>/!s/\.//g' $i 
    sed -i '/>/!s/\*//g' $i 
done


