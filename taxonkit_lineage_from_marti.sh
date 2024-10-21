#!/bin/bash
#https://bioinf.shenwei.me/taxonkit/usage/


#Name of sample, will change this with each analysis
#sample=../flush_experiment_160124/flush_experiment_marti_lca_0.1_all_levels_2024-FEB-22_10-25-17
sample=../Regular_Collections_2023/marti_assignments_lca_0.0_all_levels_2024-OCT-8_16-34-1

#Extract NCBI ID which is the 2nd field
#cut -f 2 -d $',' ${sample}.csv > ${sample}_taxaID.txt #for csv data
cut -f 2 -d $'\t' ${sample}.tsv > ${sample}_taxaID.txt #for tsv data



#Get lineages
taxonkit lineage ${sample}_taxaID.txt > ${sample}_taxaID_lineage.txt


#Fill in blanks
taxonkit reformat ${sample}_taxaID_lineage.txt -r Unassigned | cut -f 1,3 > ${sample}_taxaID_lineage_clean.txt


#Change delimiter
sed 's/;/\t/g' ${sample}_taxaID_lineage_clean.txt | awk -F'\t' 'BEGIN {OFS=","} { print $1, $2, $3, $4, $5, $6, $7, $8 }' > ${sample}_taxaID_lineage_sep.csv


#Add headers
echo "NCBI ID,kingdom,phylum,class,order,family,genus,species" > header.txt
cat header.txt ${sample}_taxaID_lineage_sep.csv > ${sample}_taxaID_lineage.csv


#Remove all the intermediate files generated
rm ${sample}_taxaID.txt ${sample}_taxaID_lineage.txt ${sample}_taxaID_lineage_clean.txt ${sample}_taxaID_lineage_sep.csv header.txt