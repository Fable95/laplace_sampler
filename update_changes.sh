#!/usr/bin/env bash

echo "Update submodule"

git submodule update --init --recursive

echo "Move changes to MP-SPDZ"

source_directory="./changes"
target_directory="./MP-SPDZ"

if [ -d "$source_directory" ]; then
    if [ -d "$target_directory" ]; then
        for item in "$source_directory"/*; do
            echo "Move '$item' to '$target_directory'"
            cp -Rf "$item" "$target_directory" 
        done
    else
        echo "'$target_directory' directory does not exist"
    fi  
else
    echo "'$source_directory' directory does not exist"
fi