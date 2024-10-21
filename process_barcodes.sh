#!/bin/bash

# This script processes Nanopore sequencing reads for multiple reference genomes and barcodes.
# It runs the calculate_minimapped_bases.py script on each barcode file and captures the output.
# The script generates a summary file for each reference genome, containing the total matching bases
# and the number of mapped reads for each barcode.

#Is set up to work when run from ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/church_farm_2023

#First draft written Oct 2024

#Set up list of reference genomes - these need to match the file paths in the genome_coverage directory
ref_genome_list=("z_tritici" "p_avenae" "p_nodorum")

# Define an array of barcode numbers using brace expansion
barcode_list=$(seq -w 01 39)

for ref_genome in "${ref_genome_list[@]}"; do

    # Output file 
    output_file="genome_coverage/${ref_genome}_coverage.tsv"
    # Create output file if it doesn't exist
    touch $output_file

    # Write the header to the output file
    echo -e "Reference Genome\tBarcode\tTotal Matching Bases\tNumber of Mapped Reads" > $output_file

    # Loop through each barcode number
    for barcode in $barcode_list; do
        # Construct the barcode file path
        barcode_file="genome_coverage/${ref_genome}/barcode${barcode}_mapped.paf"
        
        # Run the calculate_minimapped_bases.py script and capture the output
        output=$(python scripts/calculate_minimapped_bases.py "$barcode_file")
        
        # Append the output to the new file
        echo -e "$output" >> $output_file
    done

    echo "Summary written to $output_file"
done

echo "All done!"