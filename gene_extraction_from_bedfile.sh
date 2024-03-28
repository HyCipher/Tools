#/bin/bash

work_path="/BiO/Live/rooter/Downloads/synteny/bed_data"

cd $work_path

for bed_file in *.bed; do
    bedtools sort -i $bed_file
    grep 'ID=gene' $bed_file > ./gene_bed/$(basename $bed_file .bed)_gene.bed
    done
