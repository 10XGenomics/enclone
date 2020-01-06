# enclone

The mission of enclone is to:

**Find and display the clonotypes within single cell VDJ datasets:
groups of cells having the same fully rearranged common ancestor.**

enclone is part of the 10x Genomics immune profiling tools, including
Cell Ranger and Loupe, which enclone is integrated with.

To use enclone, you need to know how to run command-line tools.  This is something that 
can be learned easily, particularly if you have a friend or colleague who can help you
get started.  You do not need to be able to program, or anything of that sort.

Here are *temporary* instructions for using enclone (for internal testing purposes):

1. For now, you can run on an x86-64 linux server, or a Mac, and possibly on a Windows
box (untested).

2. You need to have the rust compiler installed.

3. Check out the branch `dj/cr-1577b` of `cellranger`, and do
```
cd cellranger/lib/rust
cargo build --release --bin enclone
```
and then put `target/release/enclone` in your path.

4. Copy the directory `enclone/test/inputs` to somewhere you can point to, or just leave it 
where it is.  These are test data you can play with.  (You can also supply your own output
from a Cell Ranger immune profiling run.)  When you read the documentation at step 6, you'll 
get to a place where you put `PRE=enclone/test/inputs` or instead with the path where you've
copied the data to.  But you need to supply `PRE` with a path that makes sense from the directory
you're working in.

5. Type `enclone help`, and read the terminal setup instructions there.

6. Type `enclone` and study the documentation shown there.

7. If you want to run the built-in tests, type
```
cargo test --release -p enclone enclone -- --nocapture
```

If you have problems, please write to us at enclone@10xgenomics.com.
