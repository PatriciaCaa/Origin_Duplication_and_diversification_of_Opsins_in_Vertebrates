### PHYLOGENTIC ANALYSIS

#!/bin/bash

# STEP 1: Reduce redundancy using CD-HIT at various similarity thresholds
# Input: AllOpsinhits.fasta (from homology identification stage)
echo "Running CD-HIT clustering at multiple thresholds..."
cd-hit -i AllOpsinhits.fasta -o AllOpsinhits_90.fasta -c 0.9
cd-hit -i AllOpsinhits.fasta -o AllOpsinhits_85.fasta -c 0.85
cd-hit -i AllOpsinhits.fasta -o AllOpsinhits_80.fasta -c 0.8
cd-hit -i AllOpsinhits.fasta -o AllOpsinhits_70.fasta -c 0.7
cd-hit -i AllOpsinhits.fasta -o AllOpsinhits_50.fasta -c 0.5

# STEP 2: Add melatonin receptors to opsin dataset 
# (Assumes you've already concatenated melatonin receptor sequences)
cat AllOpsinhits_80.fasta melatonin_receptors.fasta > AllOpsinsPlusMelatonin_80.fasta

# STEP 3: Multiple sequence alignment with MAFFT
echo "Running MAFFT alignment..."
mafft --maxiterate 1000 AllOpsinsPlusMelatonin_80.fasta > AllOpsins_80_aligned.fasta

# STEP 4: Trim poorly aligned regions using TrimAl
echo "Trimming alignment..."
trimal -in AllOpsins_80_aligned.fasta -gt 0.2 -out AllOpsins_80_aligned_trimmed.fasta

# STEP 5: Build preliminary phylogenetic tree using FastTree
echo "Constructing FastTree..."
fasttree -lg -gamma -cat 1 AllOpsins_80_aligned_trimmed.fasta > AllOpsins_fasttree.tree

# STEP 6: Model selection and final tree with IQ-TREE2 
echo "Running IQ-TREE2 with model GTR20+F+G4..."
iqtree2 -s AllOpsins_80_aligned_trimmed.fasta \
  -m GTR20+F+G4 \
  -bb 1000 \
  -alrt 1000 \
  -nt AUTO \
  -pre IQTree_Opsins_80

echo "DONE. Tree written to IQTree_Opsins_80.treefile"

