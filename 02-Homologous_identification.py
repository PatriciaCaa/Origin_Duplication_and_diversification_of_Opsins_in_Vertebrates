### HOMOLOGOUS IDENTIFICATION

#!/bin/bash

# STEP 1: Clean protein sequences in all .faa files
# Remove any stray characters (e.g., dots, asterisks)
echo "Cleaning .faa files..."
for i in *.faa; do
  sed -i '/>/!s/\.//g' "$i"
  sed -i '/>/!s/\*//g' "$i"
done

# STEP 2: Run DIAMOND BLASTp using known opsin sequences
echo "Running DIAMOND BLASTp searches..."
for i in *.faa; do
  diamond makedb --in "$i" -d "$i.db"
  diamond blastp \
    -q Opsin_sequences.fasta \
    -d "$i.db" \
    -e 1e-10 \
    --ultra-sensitive \
    -o "$i.opsin.hits" \
    -k 50 \
    -f 6
done

# STEP 3: Extract unique hit names from DIAMOND results
echo "Extracting unique sequence names from DIAMOND results..."
for i in *.opsin.hits; do
  awk '{print $2}' "$i" | sort | uniq > "$i.names"
done

# STEP 4: Extract matched sequences from original proteomes using seqkit
echo "Extracting sequences using seqkit..."
for i in *.faa; do
  seqkit grep -n -f "$i.opsin.hits.names" "$i" > "$i.Hitseq"
done

# STEP 5: Run InterProScan on extracted sequences to detect opsin domains
echo "Running InterProScan..."
for i in *.Hitseq; do
  interproscan.sh -i "$i" -f TSV -appl PRINTS,ProSitePatterns -o "$i.Interpro.tsv"
done

# STEP 6: Run Phobius to predict transmembrane domains
echo "Running Phobius predictions..."
for i in *.Hitseq; do
  perl /path/to/phobius.pl -s "$i" > "$i.PhobOut"
done

# Clean Phobius output
for i in *.PhobOut; do
  sed '1d' "$i" | awk '{print $1,$2,$3}' > "$i.fix"
done

# STEP 7: Filter for 6-8 TMD sequences and extract them
echo "Filtering sequences with 6-8 transmembrane domains..."
for i in *.PhobOut.fix; do
  grep -E ' 6 | 7 | 8 ' "$i" >> "$i.7TMDseqs"
  sed 's/\t/ /g' "$i.7TMDseqs" | awk '{print $1}' > "$i.7TMDseqs.names"
done

for i in *.faa.Hitseq; do
  base="${i%.faa.Hitseq}"
  seqkit grep -n -f "$base.faa.Hitseq.PhobOut.fix.7TMDseqs.names" "$i" > "$i.Opsin7TMD.fasta"
done

# STEP 8: Merge all extracted 7TMD sequences into a single file
echo "Merging all sequences into AllOpsinhits.fasta..."
cat *.Opsin7TMD.fasta > AllOpsinhits.fasta

echo "DONE"





