#!/bin/bash

# This is the installation and update script for enclone.  For instructions on how to
# run it, please see bit.ly/enclone.  A few more details are here.

# This reuses code from the installation script for the rust language.
# [TO ACKNOWLEDGE ELSEWHERE.]

# This script expects a single argument, which is small, medium or large, depending on how much
# data is to be downloaded.  
#
# If you run it a second time and forget the size, it will use the same size as last time,
# and state that.
#
# The script assumes that curl is installed on your computer.  Some linux computers may not
# have this.

first=$1

#  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

main() {

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 1. Test for existence of needed system commands.
    #
    #    [SHOULD WE MAKE THIS SCRIPT WORK WITH EITHER CURL OR WGET?]

    need_cmd uname
    need_cmd curl
    need_cmd mkdir
    need_cmd chmod
    need_cmd awk

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 2. Determine if this is a Linux or Mac box; fail if it is not one of these two.

    local _ostype
    _ostype="$(uname -s)"
    if [ "$_ostype" != Linux ] && [ "$_cputype" != Darwin ]; then
        echo
        echo "enclone install script fails because operating system type ${_ostype}" \
            "is unknown."
        echo "If you're stuck please ask for help by emailing enclone@10xgenomics.com."
        echo
        exit 1
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 3. Test input arguments.

    if [ "$first" != small ] && [ "$first" != medium ] && [ "$first" != large ]; then
        echo
        echo "enclone install script fails because the first argument to it is missing" \
            "or unrecognized."
        echo "The only allowed arguments are small, medium and large."
        echo "If you're stuck please ask for help by emailing enclone@10xgenomics.com."
        echo
        exit 1
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 4. Determine if the local enclone executable is current.

    local _current_version _enclone_is_current
    _enclone_is_current=false
    if test -f "$HOME/bin/enclone"; then
        local _local_version
        if test -f "$HOME/enclone/version"; then
            _local_version=$(cat $HOME/enclone/version)
            _current_version=$(curl -sI \
                https://github.com/10XGenomics/enclone/releases/latest/download/enclone_linux | \
                grep "^location:" | tr '/' ' ' | cut -d ' ' -f9)
            if [ "$_local_version" == "$_current_version" ]; then
                echo "The local version of enclone is current so not downloading executable."
                _enclone_is_current=true
            fi
        fi
    fi

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 5. Make directory ~/bin if needed and download the appropriate enclone executable into it.

    if [ "$_enclone_is_current" = false ]; then
        mkdir -p ~/bin
        cd ~/bin
        if [ "$_ostype" = Linux ]; then
            echo "Downloading the Linux version of the latest enclone executable."
            curl -L https://github.com/10XGenomics/enclone/releases/latest/download/enclone_linux \
                --output enclone
        fi
        if [ "$_ostype" = Darwin ]; then
            echo "Downloading the Mac version of the latest enclone executable."
            curl -L https://github.com/10XGenomics/enclone/releases/latest/download/enclone_macos \
                --output enclone
        fi
        # set execute permission on the enclone executable
        chmod +x enclone
        # record local version
        mkdir -p ~/enclone
        echo "$_current_version" > ~/enclone/version
    fi

    exit

    #  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    # 6. Add ~/bin to path if needed.
    #
    #    This is complicated because some versions of Linux use the file .bash_profile,
    #    and some use .profile.
    #    If the instructions here don't work, this post may be helpful:
    #    https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path.

    if [[ ":$PATH:" != *":$HOME/bin:"* ]]; then
        test -r ~/.bash_profile && echo 'PATH=~/bin:$PATH' >> ~/.bash_profile || \
            echo 'PATH=~/bin:$PATH' >> ~/.profile
    fi

}

#  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

need_cmd() {
    if ! check_cmd "$1"; then
        echo
        echo "enclone install script fails because the command $1 was not found."
        echo "If you're stuck please ask for help by emailing enclone@10xgenomics.com."
        echo
        exit 1
    fi
}

check_cmd() {
    command -v "$1" > /dev/null 2>&1
}

main
