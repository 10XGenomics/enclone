# enclone

`enclone` is a computational tool for studying immunology.  The mission of `enclone` is to:

**Find and display the clonotypes within single cell VDJ datasets:
groups of T and B cells having the same fully rearranged common ancestor.**

The following diagram shows what a _clonotype_ is for B cells.  The same applies for T cells,
but things are simpler because T cells do not have somatic hypermutation (SHM).

<img src="img/what_is_a_clonotype.png" alt="what is a clonotype" title="what is a clonotype" />

Each cell in a clonotype is typically represented by two or three chains.  Such information can
only be obtained from _single cell_ data!  From such data, clonotypes can be computationally
approximated, with high accuracy.  The method we use for this is describe briefly in the online
documentation for enclone, and will be described separately in more detail.

___________________________________________________________________________________________________


`enclone` is beta software.  Binary executables from linux, mac and windows can be directly
downloaded from this page, as can sample 10x Genomics datasets.  And then you're off and
running!  To use `enclone`, you need to know how to run command-line tools.  This is something that 
can be learned easily, particularly if you have a friend or colleague who can help you
get started.  You do not need to be able to program, or anything of that sort.

`enclone` is part of the [10x Genomics](https://www.10xgenomics.com/) immune 
profiling toolkit, including
[Cell Ranger and Loupe](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome), 
with which `enclone` will be integrated (later).

___________________________________________________________________________________________________


The example below shows how `enclone` prints out clonotypes.  This is something you will need
to study in order to use enclone successfully.  `enclone` comes with extensive online 
documentation, and because you can easily play with the sample datasets, and gradually
figure out how it all works.

<p><img src="img/enclone_annotated_example.svg" alt="enclone annotated example" title="enclone annotated example" /></p>

By default, `enclone` prints clonotypes in this human-readable form.  You can also instruct
`enclone` to print clonotypes in machine-readable forms that are suitable for input to other
programs.

___________________________________________________________________________________________________


Here are *temporary* instructions for using `enclone` (for internal testing purposes):

1. For now, you can run on an x86-64 linux server, or a Mac, and possibly on a Windows
box (untested).

2. You need to have the Rust compiler installed. Detailed instructions on how to do this
can be found [here](https://www.rust-lang.org/tools/install). You can confirm that you 
have successfully installed the Rust compiler by running `rustc --version`.

3. Clone the `enclone` repository and build `enclone` using Cargo (which comes with Rust) by running:
```
git clone git@github.com:10XDev/enclone.git
cargo build --release --bin enclone
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; and then put `target/release/enclone` in your path. This can be done in multiple ways, including running
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`PATH=$PATH:~/where/you/cloned/enclone/target/release/enclone`. 
See also [this helpful post on StackExchange](https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Compilation takes 8-10 minutes on a 2017 MacBook Pro with a dual-core i7 and 5-7 minutes on a similar Linux machine. 

4. Copy the directory `enclone/test/inputs` to somewhere you can point to, or just leave it 
where it is.  These are test data you can play with; you can also supply your own output
from a Cell Ranger immune profiling run (so long as there is an `all_contig_annotations.json` output). 
When you read the documentation at step 6, you'll get to a place where you put `PRE=enclone/test/inputs` 
or instead with the path where your copied data reside.  But you need to supply `PRE` with a path that 
makes sense relative to your working directory.

5. Type `enclone help`, and read the terminal setup instructions there.

6. Type `enclone` and study the documentation shown there.

7. If you want to run the built-in tests, type
```
cargo test --release -- --nocapture
```

If you have problems, please write to us at enclone@10xgenomics.com.

