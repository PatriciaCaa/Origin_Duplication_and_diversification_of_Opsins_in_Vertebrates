
### BAsh Script to download genome datasets in a automated mode using pre-selected accession numbers

#!/bin/bash

while read line

do
    datasets download genome accession ${line} --include protein #ensuring the sequences are included
    unzip ncbi_dataset.zip
    rm README.md

# Move the extracted protein sequence file (protein.faa) to a new location and rename the file to include the accession number
    mv ncbi_dataset/data/${line}/protein.faa ${line}_protein.faa
    rm -rf ncbi_dataset* ## Remove the extracted dataset directory
done < acc.txt



