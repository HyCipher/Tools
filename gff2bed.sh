#!/bin/bash

# Specify the directory the file in
file_path="/BiO/Live/rooter/Downloads/synteny/bed_data"

# Change the gff file to bed file
for file in "$file_path"/*.gff; do
    gff2bed < $file > "$(basename $file .gff).bed"
    done