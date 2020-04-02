# Fetching test datasets for enclone

Here are two methods that should do the same thing.  First `cd` so that you're in your
home directory.

1.  `svn export https://github.com/10XGenomics/enclone/trunk enclone`

2.  `curl -L https://github.com/10XGenomics/enclone/archive/master.tar.gz | tar zx enclone-master/test; mv enclone_master enclone`

If, perchance, neither `svn` nor `curl` are installed on your computer, you will need to 
install them.  We will install instructions for this later.
