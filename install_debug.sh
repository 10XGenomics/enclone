#!/bin/bash

printf "\nSome diagnostic information will be printed out below.\n"

printf "\n1. Determining which shell you are using: $SHELL.\n"

printf "\n2. Show the permissions on ~/bin/enclone:\n\n"
ls -l ~/bin/enclone

printf "\n3. Testing for existence of various initialization files in your home directory\n"
printf "   and for each such file, if present, whether it sets your path.\n\n"
cd
NEWLINE=1
    
# Decide which initialization files to test.
    
if [ "$SHELL" == "/bin/tcsh" ]; then
    files=( ".tcshrc" ".login" ".cshrc" )
elif [ "$SHELL" == "/bin/zsh" ]; then
    files=( ".zshenv" ".zprofile" ".zshrc" ".zlogin" )
elif [ "$SHELL" == "/bin/bash" ]; then
    files=( ".profile" ".bash_login" ".bash_profile" ".bashrc" )
elif [ "$SHELL" == "/bin/sh" ]; then
    files=( ".profile" ".bash_login" ".bash_profile" ".bashrc" )
else
    files=( ".profile" \
        ".zshenv" ".zprofile" ".zshrc" ".zlogin" \
        ".bash_login" ".bash_profile" ".bashrc" \
        ".tcshrc" ".login" ".cshrc" )
fi
    
# Test initialization files.
    
for i in "${files[@]}"
do
    ls -l $i >& /dev/null
    if [ "$?" -eq "0" ]; then
         if [ "$NEWLINE" -eq "0" ]; then
             echo ""
             NEWLINE=1
         fi
         printf "$i: present\ntesting it for setting path\n"
         cat $i | grep -i PATH | grep -v "^#" | uniq
         echo ""
    else
         echo "$i: absent"
         NEWLINE=0
    fi
done
if [ "$NEWLINE" -eq "0" ]; then
    echo ""
