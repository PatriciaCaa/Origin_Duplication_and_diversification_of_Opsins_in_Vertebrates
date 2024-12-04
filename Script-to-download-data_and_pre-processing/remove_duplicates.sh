### Script used to remove duplicated files


#!/bin/bash

#Path to the directory containing my .faa files
directory="/home/vpc5/IRP/singlecell_genome_data/Singlec_genomes_copy"

#Iterate through all .faa files in the directory
for file in "$directory"/*.faa; do
    if [ -f "$file" ]; then
        #output file name
        output_file="${file%.faa}_deduplicated.faa"

        #seqkit rmdup command - this will remove any duplicated files
        echo "Processing $file -> $output_file"
        seqkit rmdup -n "$file" -o "$output_file"

        if [ $? -eq 0 ]; then
            echo "Processed $file -> $output_file"
        else
            echo "Failed to process $file"
        fi
    else
        echo "No .faa files found in the directory."
    fi
done
