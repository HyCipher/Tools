#!/bin/bash

# -------------------------------------------------------------------------------
# Customize section
# Get the number of file lines
file="/BiO/Live/rooter/Downloads/synteny/ortho_method/proteomes/test.txt"
output_file="/BiO/Live/rooter/Downloads/synteny/tools/1.txt"
# -------------------------------------------------------------------------------

while IFS= read -r line; do
    # -------------------------------------------------------------------------------
    # Customize section
    count=$(wc -l "/BiO/Live/rooter/Downloads/synteny/ortho_method/proteomes/${line}_protein/primary_transcripts/OrthoFinder/Results_Mar20/Orthogroups/Orthogroups_UnassignedGenes.tsv" | awk '{print $1}')
    # -------------------------------------------------------------------------------
    echo "$line,$count" >> "$output_file"
done < "$file"
