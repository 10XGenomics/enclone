# Fetching test datasets for enclone

Here are two methods that should do the same thing.  First `cd` so that you're in your
home directory.

1.  `svn export https://github.com/10XGenomics/enclone/trunk enclone`

2.  `curl -L https://github.com/10XGenomics/enclone/archive/master.tar.gz | tar zx enclone-master/test; mv enclone_master enclone`

3.  `wget -O- https://github.com/10XGenomics/enclone/archive/master.tar.gz | tar zx enclone-master; mv enclone-master enclone`

If, perchance, neither `svn` nor `curl` nor `wget` are installed on your computer, you will need 
to install one of them.  We will add installation instructions for this later.

Or maybe this will help:

<img align="left" src="https://imgs.xkcd.com/comics/universal_install_script.png" alt="universal install script" title="universal install script" />
