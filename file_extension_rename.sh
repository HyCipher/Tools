#!/bin/bash
# -------------------------------------------------------------------------------
# Customize section
# Specify the file extension to be modified指定要修改的文件名结尾
old_extension=".fna"
new_extension=".fa"

# Specify the directory the file in
target_directory="/BiO/Live/rooter/Downloads/synteny/ortho_method/fasta_data"
# -------------------------------------------------------------------------------

# Change into the target directory
cd "$target_directory" || exit
# Iterate through all files with old extensions in the target directory
for file in *"$old_extension"; do
    # Check if the file is an ordinary file
    if [ -f "$file" ]; then
        # Build a new name for file
        new_name="${file%"$old_extension"}$new_extension"

        # Rename
        mv "$file" "$new_name"

        # Report
        echo "Successul rename '$file' to '$new_name'"
    fi
done

echo "All file extension modifications successfully in batches"
