#!/bin/bash

# Parsing command-line options. 
# Show help message if user calls 'rename_files -h'
# If any other option is provided, give error message.
# -->   The colon before the list of allowed options is needed to 
#       disable default error handling (which you should) and use '\?'
while getopts ":h" opt; do
    case ${opt} in
        h ) echo "Usage: rename_files [-h] [oldname] [newname]"; exit;
            ;;
        \? ) echo "Nope. Invalid option provided."; exit;
            ;;
        esac
done

# "The [[ ]] construct is more versatile Bash version of [].
# This is the 'extended test command'. "
if [ $# -eq 2 ]; then
    OLDNAME="$1"
    NEWNAME="$2"
else 
    OLDNAME="TowerDumper"
    NEWNAME="MyCaloEvaluator"
fi

# __________ 1. Global search and replace within files. _________
# find [path...] [expression]
# Test:     '-type f': File is of type 'regular file'.
# Action:   '-exec': All following args are considered in command until ';'
#           --> The string '{}' is replaced by the current file name being 
#               processed everywhere it occurs in the args to the command. 
#           --> 'exec COMMAND {} +': Runs command on the selected files, but the 
#               command line is build by appending each selected file name at the end.
echo "Executing global search and replace: ${OLDNAME} --> ${NEWNAME}"
find . ! -name '*rename_files.sh*' -type f -exec sed -i "s/${OLDNAME}/${NEWNAME}/Ig" {} +
if grep -iq --exclude="rename_files.sh" "${OLDNAME}" **; then
    echo "Error: Found occurrence of ${OLDNAME} after search and replace. Terminating."
    exit
fi

# Only for files with the OLDNAME prefix. 
for file in *${OLDNAME}*; do
    if [ file == "rename_files.sh" ]; then 
        echo "woops lulz"
        exit
    fi
    
    # Syntax: ${string##substr}
    # -->   DELETES longest match of substr from FRONT of string; returns what is leftover. 
    # -->   Here, extension is set to whatever comes after the first period in a filename, 
    #       since everything before, and including, the period is deleted.
    # NOTE: Not actually used but leaving for future reference. 
    extension="${file##*.}"

    # Syntax: ${string//substr/replacement}
    # -->   Replace all matches of $substr with $replacement.
    echo "Renaming $file to ${file//${OLDNAME}/${NEWNAME}}"
    mv "${file}" "${file//${OLDNAME}/${NEWNAME}}"
done
