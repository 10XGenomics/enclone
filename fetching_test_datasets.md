# Fetching test datasets for enclone

Here are three methods that should do the same thing.  First `cd` so that you're in your
home directory.  Then do <b>one</b> of the following:

1.  `svn export https://github.com/10XGenomics/enclone/trunk enclone`

2.  `curl -L https://github.com/10XGenomics/enclone/archive/master.tar.gz | tar zx enclone-master/test; mv enclone_master enclone`

3.  `wget -O- https://github.com/10XGenomics/enclone/archive/master.tar.gz | tar zx enclone-master; mv enclone-master enclone`

# Installing curl
If, perchance, neither `svn` nor `curl` nor `wget` are installed on your computer, you will need 
to install one of them.  There are a few ways to do this depending on which flavor of Unix you use. Here are some ways to get `curl` onto your machine.

### Debian/Ubuntu
For these machines, use <b>one</b> of the following:

1. `apt-get install curl`

2. `sudo apt-get install curl`

3. `sudo apt install curl`

### Mac OS X
To get the latest `curl` on OS X, you will likely want to use the latest [Homebrew](https://brew.sh):

1. `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"`
2. `brew install curl`

### CentOS/Fedora/RHEL
For these machines, use the following:

`yum install curl`

### Other flavors
If you're using another flavor of Unix, these may be helpful:

1. <b>OpenSUSE:</b> `zypper install curl`

2. <b>ArchLinux:</b> `pacman -Sy curl`

# Installing wget
Have something against `curl` but don't know how to install `wget`? This section is for you.

### Debian/Ubuntu
For these machines, use <b>one</b> of the following:

1. `apt-get install wget`

2. `sudo apt-get install wget`

3. `sudo apt install wget`

### Mac OS X
To get the latest `wget` on OS X, you will likely want to use the latest [Homebrew](https://brew.sh):

1. `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"`
2. `brew install wget`

### CentOS/Fedora/RHEL
For these machines, use the following:

`yum install wget`

### Other flavors
If you're using another flavor of Unix, these may be helpful:

1. <b>OpenSUSE:</b> `zypper install wget`

2. <b>ArchLinux:</b> `pacman -Sy wget`

# Installing svn
If you'd rather use `svn` instead of `wget` or `curl`, feel free! If you don't know how to install `svn`, this section
should help you do so.

### Ubuntu
For these machines, run the following commands <b>in order</b>:

1. `apt-get install subversion`

2. `apt-get install libapache2-svn`

### Debian
For these machines, run the following commands <b>in order</b>:

1. `apt-get install subversion`

2. `apt-get install libapache2-mod-svn`

### Mac OS X
To get the latest `svn` on OS X, you will likely want to use the latest [Homebrew](https://brew.sh):

1. `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"`
2. `brew install subversion`

### Fedora
For Fedora machines, run `yum install subversion`.

### CentOS/RHEL
For these machines, run the following commands <b>in order</b>:

1. `yum install subversion`

2. `yum install mod_dav_svn`

### Other flavors
If you're using another flavor of Unix, these may be helpful:

1. <b>OpenSUSE:</b> `zypper install subversion` and then `zypper install subversion-server`

2. <b>ArchLinux:</b> `pacman -Sy subversion`

# I still can't install curl, wget, or svn!
If all else fails, maybe these instructions from xkcd will help:

<img align="left" src="https://imgs.xkcd.com/comics/universal_install_script.png" alt="universal install script" title="universal install script" />

Well seriously, if you're really stuck, please write to us at enclone@10xgenomics.com.
