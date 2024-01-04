#!/bin/bash

source_folder="/broad/IDP-Dx_work/nirmalya/local/bin"
target_folder="/home/unix/hsethura/RNA-Seq-Pipeline/broad/IDP-Dx_work/nirmalya/local/bin"

for file in "$source_folder"/*; do
    if [ -f "$file" ]; then
            ln -s "$file" "$target_folder/$(basename "$file")"
    elif [ -d "$file" ]; then
            ln -s "$file" "$target_folder/$(basename "$file")"
    fi
done

