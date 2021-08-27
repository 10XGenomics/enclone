#!/bin/bash

# Note that we call this script using "bash -c", and that's what determines the shell,
# not the first line of this file.

# This is the installation and update script for enclone.  For instructions on how to
# run it, please see bit.ly/enclone.  A few more details are here.
#
# Run shellcheck if you change this script: you'll need to sort wheat from chaff.

#  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

# This reuses code from the installation script for the rust language.
#
# This script expects a single argument, which is small, medium or large, depending on how much
# data is to be downloaded.  
#
# If you run it a second time and forget the size, it will use the same size as last time,
# and state that.
#
# Note that version15 is hardcoded!

#  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

# Get command line arguments.  The second and third arguments are for testing.
# The third argument can be force_wget.

size=$1
if ! [ -z "$2" ]; then
    HOME=$2
fi
if ! [ -z "$3" ]; then
    MODE=$3
fi

#  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

main() {

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 1. Set up; test for existence of needed system commands.  
    #
    #    We require only one of curl or wget.  The reason for not requiring curl is that at
    #    the time of writing this script, the standard Ubuntu install did not include curl,
    #    and so it is possible that someone would not have curl.
    #
    #    We do not use svn, because it is no longer available by default on MacOS.
    #
    #    Sadly, no longer is git.  We therefore include a special test for it.

    need_cmd date
    STARTTIME=$(date +%s)
    # special test for git
    git --version >& /dev/null
    if ! [ "$?" -eq "0" ]; then
        printf "\nIt would appear that you do not have the command line tool git installed.\n"
        printf "This is a common problem.  To solve it, please type\n"
        printf "xcode-select --install\n"
        printf "and follow the instructions, and then rerun the enclone installation command.\n\n"
        exit 1
    fi
    # We used to have "set -e" here to force exit upon error.  However that is inconsistent with
    # several parts of the code that allow for a command to fail and then check status.
    # There has to have been a reason why we added set -e.  If we figure out what it was, we can
    # protect the specific parts of code that need protection.
    need_cmd awk
    need_cmd chmod
    need_cmd grep
    need_cmd mkdir
    need_cmd uname
    need_cmd zcat
    local _have_curl
    _have_curl=false
    if check_cmd curl; then
        if [ "$MODE" != force_wget ]; then
            _have_curl=true
        fi
    fi
    if ! $_have_curl && ! check_cmd wget; then
        printf "\nenclone installation failed because neither the command curl nor the\n"
        printf "command wget could be found.  This is strange and unexpected.\n"
        printf "If you're stuck please ask for help by emailing enclone@10xgenomics.com.\n\n"
        exit 1
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 2. Determine if this is a Linux or Mac box; fail if it is not one of these two.

    local _ostype
    _ostype="$(uname -s)"
    if [ "$_ostype" != Linux ] && [ "$_ostype" != Darwin ]; then
        echo
        echo "enclone install script fails because operating system type ${_ostype}" \
            "is unknown."
        echo "If you're stuck please ask for help by emailing enclone@10xgenomics.com."
        echo
        exit 1
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 3. Get requested size.

    if [ "$size" != small ] && [ "$size" != medium ] && [ "$size" != large ]; then
        printf "\nTo install or update enclone, please supply the single argument SIZE to the\n"
        printf "curl command shown on bit.ly/enclone.  The argument SIZE can be small, medium "
        printf "or large.\n"
        echo "If you're stuck please ask for help by emailing enclone@10xgenomics.com."
        echo
        exit 1
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 4. Determine if datasets are current.
    #
    #    Because there has been only one release of the large dataset collection, if it
    #    was downloaded, then it is current.

    local _datasets_small_current _datasets_medium_current _datasets_large_current
    local _datasets_small_checksum_master _datasets_small_checksum_local
    local _datasets_medium_checksum_master _datasets_medium_checksum_local
    _datasets_small_current=false
    _datasets_medium_current=false
    _datasets_large_current=false
    raw_repo=https://raw.githubusercontent.com/10XGenomics/enclone
    if [ "$size" = small ]; then
        if [ $_have_curl ] && [ "$MODE" != force_wget ]; then
            _datasets_small_checksum_master=$(curl -s $raw_repo/master/datasets_small_checksum)
        else
            _datasets_small_checksum_master=$(wget -q $raw_repo/master/datasets_small_checksum -O -)
            if ! [ "$?" -eq "0" ]; then
                printf "\nfailed: wget -q $raw_repo/master/datasets_small_checksum\n"
                printf "This is strange and unexpected.\n"
                echo "If you're stuck please ask for help by emailing enclone@10xgenomics.com."
                echo
                exit 1
            fi
        fi
    fi
    if test -f "$HOME/enclone/datasets_small_checksum"; then
        _datasets_small_checksum_local=$(cat $HOME/enclone/datasets_small_checksum)
        if [ "$_datasets_small_checksum_local" = "$_datasets_small_checksum_master" ]; then
            _datasets_small_current=true
        fi
    fi
    if [ "$size" = medium ] || [ "$size" = large ]; then
        raw_master=$raw_repo/master
        if $_have_curl; then
            _datasets_medium_checksum_master=$(curl -s $raw_master/datasets_medium_checksum)
        else
            _datasets_medium_checksum_master=$(wget -q $raw_master/datasets_medium_checksum -O -)
            if ! [ "$?" -eq "0" ]; then
                printf "\nfailed: wget -q $raw_repo/master/datasets_medium_checksum\n"
                printf "This is strange and unexpected.\n"
                echo "If you're stuck please ask for help by emailing enclone@10xgenomics.com."
                echo
                exit 1
            fi
        fi
    fi
    if test -f "$HOME/enclone/datasets_medium_checksum"; then
        _datasets_medium_checksum_local=$(cat $HOME/enclone/datasets_medium_checksum)
        if [ "$_datasets_medium_checksum_local" = "$_datasets_medium_checksum_master" ]; then
            _datasets_medium_current=true
        fi
    fi
    if test -f "$HOME/enclone/datasets2/download_complete"; then
        _datasets_large_current=true
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 5. Determine if the local enclone executable is current.  
    #
    #    This is quite hideous.  There must be a better way.

    local _current_version _enclone_is_current _is_update
    repo=https://github.com/10XGenomics/enclone
    if $_have_curl; then
        _current_version=$(curl -sI $repo/releases/latest/download/enclone_linux | \
            grep -i "^location:" | tr '/' ' ' | cut -d ' ' -f9)
    else
        _current_version=$(wget --server-response --max-redirect=0 \
            $repo/releases/latest/download/enclone_linux |& \
            grep -i " location:" | tr '/' ' ' | cut -d ' ' -f11)
    fi
    _enclone_is_current=false
    if test -f "$HOME/bin/enclone"; then
        _is_update=true
        local _local_version
        if test -f "$HOME/enclone/version"; then
            _local_version=$(cat $HOME/enclone/version)
            if [ "$_local_version" == "$_current_version" ]; then
                printf "\nThe local version of enclone is current so not downloading executable.  "
                printf "Both versions are $_local_version.\n"
                _enclone_is_current=true
            fi
        fi
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 6. Make directory ~/bin if needed and download the appropriate enclone executable into it.

    cd $HOME
    mkdir -p bin
    mkdir -p enclone
    if [ "$_enclone_is_current" = false ]; then
        cd bin
        if [ "$_ostype" = Linux ]; then
            printf "\nDownloading the Linux version of the latest enclone executable.\n\n"
            if $_have_curl; then
                curl -s -L $repo/releases/latest/download/enclone_linux --output enclone
            else
                wget -q $repo/releases/latest/download/enclone_linux -O enclone
            fi
        fi
        if [ "$_ostype" = Darwin ]; then
            printf "\nDownloading the Mac version of the latest enclone executable.\n\n"
            if $_have_curl; then
                curl -s -L $repo/releases/latest/download/enclone_macos --output enclone
            else
                wget -q $repo/releases/latest/download/enclone_macos -O enclone
            fi
        fi
        echo "Done downloading the enclone executable."
        # set execute permission on the enclone executable
        chmod +x enclone
        cd ..
        # record local version
        echo "$_current_version" > enclone/version
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 7. Add ~/bin to bash path if needed.
    #
    #    This does nothing if you already have ~/bin in your path.
    #
    #    This is complicated because some versions of Linux use the file .bash_profile,
    #    and some use .profile.
    #    If the instructions here don't work, this post may be helpful:
    #    https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path.
    #
    #    Note for this section and the next: putting ~/bin in the user's path will not necessarily
    #    work.  Instead they should have <home dir>/bin, as in the commands below.

    if [[ ":$PATH:" != *":$HOME/bin:"* ]]; then
        test -r .bash_profile && echo 'PATH=$HOME/bin:$PATH' >> .bash_profile || \
            echo 'PATH=$HOME/bin:$PATH' >> .profile
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 8. Add $HOME/bin to zsh path if needed.
    #
    #    (a) If .zshrc exists and we have not already added $HOME/bin to the PATH in it, do so.
    #    (b) If .zshrc does not exist but the user shell is zsh, add $HOME/bin as above.

    if [ -f .zshrc ]; then
        if [[ `cat .zshrc` != *"export PATH=$HOME/bin:"* ]]; then
            echo 'export PATH="$HOME/bin:$PATH"' >> .zshrc
        fi
    elif [[ "$SHELL" == "/bin/zsh" ]]; then
        echo 'export PATH="$HOME/bin:$PATH"' > .zshrc
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 9. Record size.

    echo "$size" > enclone/size

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 10. Download small data.  

    if [ "$size" = small ]; then
        if [ "$_datasets_small_current" = false ]; then
            printf "\nDownloading small version of datasets.\n"
            printf "This seems to take roughly five seconds, even over home wireless,\n"
            printf "however, you might have a slower connection.\n\n"
            mkdir -p enclone/datasets
            rm -rf enclone/datasets/123085
            cd enclone/datasets
            mkdir -p 123085/outs
            cd 123085/outs
            json="all_contig_annotations.json.lz4"
            url="https://github.com/10XGenomics/enclone-data/raw/master/big_inputs/version15/123085/outs/$json"
            if $_have_curl; then
                curl -s -L $url -O
            else
                wget -q $url
            fi
            cd ../../../..
            echo "$_datasets_small_checksum_master" > enclone/datasets_small_checksum
            printf "Done with that download.\n"
        else
            printf "\nSmall version of datasets already current so not downloading.\n"
        fi
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 11. Download medium data.  
    #
    #     This is not optimal, because if anything changed, all the files get re-downloaded.

    if [ "$size" = medium ] || [ "$size" = large ]; then
        if [ "$_datasets_medium_current" = false ]; then
            echo
            if [ "$size" = medium ]; then
                echo "Downloading medium version of datasets."
            fi
            if [ "$size" = large ]; then
                echo "Downloading medium version of datasets (as part of large)."
            fi
            printf "This seems to take roughly thirty seconds, even over home wireless,\n"
            printf "however, you might have a slower connection.\n\n"
            rm -rf enclone/datasets enclone/version15
            cd enclone
            git clone --depth=1 https://github.com/10XGenomics/enclone-data.git
            mv enclone-data/big_inputs/version15 datasets
            rm -rf enclone-data
            cd ..
            echo "$_datasets_medium_checksum_master" > enclone/datasets_medium_checksum
            printf "Done with that download.\n"
            # Remove a funny-looking directory, which is used by enclone only to test if 
            # weird unicode characters in a path will break it.
            rm -rf enclone/datasets/█≈ΠΠΠ≈█
        else
            printf "\nMedium version of datasets already current so not downloading them.\n"
        fi
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 12. Download large data.  

    if [ "$size" = large ]; then
        if [ "$_datasets_large_current" = false ]; then
            printf "\nDownloading large version of datasets.\n"
            printf "This seems to take roughly one to three minutes, even over home wireless,\n"
            printf "however, you might have a slower connection.\n\n"
            cd enclone
            rm -rf datasets2
            aws=https://s3-us-west-2.amazonaws.com
            if $_have_curl; then
                curl -s $aws/10x.files/supp/cell-vdj/enclone_data_2.0.tar.gz -O
            else
                wget -q $aws/10x.files/supp/cell-vdj/enclone_data_2.0.tar.gz
            fi
            cat enclone_data_2.0.tar.gz | zcat | tar xf -
            rm enclone_data_2.0.tar.gz
            mv enclone_data_2.0 datasets2
            cd ..
            touch enclone/datasets2/download_complete
            printf "Done with that download.\n"
        else
            printf "\nLarge version of datasets already current so not downloading them.\n"
        fi
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 13. Test to see if the user would get enclone, and the right version, from the
    #     command line.

    printf "\ntesting for availability of enclone by asking enclone to print its version\n\n"
    #
    # 1. Debugging is started if `$SHELL -c "enclone --version"` fails or returns the wrong version.
    # 2. Print the shell that is being used.
    # 3. Print the path that is being used by the given shell.
    # 4. Define a list of initialization files, as a function of the shell.
    # 5. For each initialization file, test if it exists, and if it exists, check for path setting.
    #
    ok=0
    $SHELL -c "enclone --version" 
    if [ "$?" -eq "0" ]; then
        available_version=$($SHELL -c "enclone --version")
        available_version=v$(echo $available_version | tr ' ' '\n' | head -1)
        if [ "$_current_version" == "$available_version" ]; then
            printf "\nGood, the version of enclone that you would get from the command line is\n"
            printf "the version $_current_version that was just downloaded.\n"
            ok=1
        else
            printf "\nIt would appear that enclone is available as a command, as determined by\n"
            printf "your current shell and shell initialization files, but your configuration\n"
            printf "does not give you the enclone version $_current_version that was just "
            printf "downloaded.\n"
        fi
    else
        printf "\nIt would appear that enclone is not available as a command, as determined by\n"
        printf "your current shell and shell initialization files.\n"
    fi
    if [ "$ok" -eq "0" ]; then
        printf "\nSome diagnostic information will be printed out below.\n"
        printf "\n1. Determining which shell you are using: $SHELL.\n"
        printf "\n2. Determining the path defined by your shell.\n\n"
        $SHELL -c "echo $PATH"
        printf "\n3. Show the permissions on ~/bin/enclone:\n\n"
        ls -l ~/bin/enclone
        printf "\n4. Attempt to execute ~/bin/enclone directly:\n\n"
        $SHELL -c "$HOME/bin/enclone --version"
        printf "\n5. Show the output of which enclone.\n\n"
        which enclone
        printf "\n6. Testing for existence of various initialization files in your home directory\n"
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
        fi
        printf "🌹 As indicated above, something has gone awry with your enclone installation. 🌹\n"
        printf "🌹                                                                             🌹\n"
        printf "🌹 Please cut and paste what is in your terminal window, as text, starting with🌹\n"
        printf "🌹 the curl command that you typed to install enclone, and send it to          🌹\n"
        printf "🌹 enclone@10xgenomics.com.  We will try to diagnose the problem!              🌹\n"
        echo
        exit 1
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 14. Succesfully completed.

    ENDTIME=$(date +%s)
    echo
    if [ "$_is_update" = true ]; then
        echo "enclone update took $(($ENDTIME - $STARTTIME)) seconds."
    else
        echo "enclone installation took $(($ENDTIME - $STARTTIME)) seconds."
    fi
    printf "\n"
    printf "🌸 If you CLOSE this terminal window and open a new one, then enclone will be     🌸\n"
    printf "🌸 in your executable path.  Otherwise enclone may not be found when you type it. 🌸\n"
    printf "\nAll done, have a lovely day!\n\n"

}

#  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

need_cmd() {
    if ! check_cmd "$1"; then
        printf "\nenclone installation faileds because the command $1 was not found.\n"
        printf "If you're stuck please ask for help by emailing enclone@10xgenomics.com.\n"
        printf "It is possible that we can rewrite the script to not use $1.\n\n"
        exit 1
    fi
}

check_cmd() {
    command -v "$1" > /dev/null 2>&1
}

main
