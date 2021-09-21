#!/bin/bash

# This prints some debugging information that may be helpful for debugging enclone installation
# problems.
#
# At one time we had this as part of code in install.sh that attempted to detect installation
# failures.  This did not work reliably.  Pretty much, there were two lessons of that attempt:
#
# 1. It is not possible for a shell script to determine what would happen if a command was typed
#    from a fresh terminal window.  You can get close, but not all the way there.
#
# 2. Every line in a shell script should be treated as a liability and dangerous.  Do not add
#    lines unless they are absolutely needed.

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
