# enclone

`enclone` is a computational tool for studying the adaptive immune system in humans and
other vertebrate species.  

It uses data from [10x Genomics](https://www.10xgenomics.com/) which permit the precise
characterization of individual immune cells.

For the general audience: every time you get sick, your body mounts an immune response involving
selective amplification of immune cells and mutations within them.  From a biological sample 
(such as a blood draw), enclone allows you to see the history of individual immune cells within a 
biological sample (such as a blood draw or biopsy): how they evolved in response to antigens, 
including those in viruses, bacteria and tumors.  Understanding this is a highly active 
research area with many mysteries that will take years if not decades to fully sort out, with 
profound implications for biology and medicine.  Welcome to the frontier!

`enclone` and this page are designed for immunologists, but if
you're simply curious, there's nothing to stop you from downloading and playing with it.

The mission of `enclone` is to:

**Find and display the clonotypes within single cell VDJ datasets:
groups of T and B cells having the same fully rearranged common ancestor.**

The following diagram shows what a _clonotype_ is for B cells.  The same applies for T cells,
but things are simpler because T cells do not have somatic hypermutation.

<img src="img/what_is_a_clonotype.png" alt="what is a clonotype" title="what is a clonotype" />

Each cell in a clonotype is typically represented by two or three chains.  Such information can
only be obtained from _single cell_ data!  From such data, clonotypes can be computationally
approximated, with high accuracy.  The method we use for this is described briefly in the online
documentation for `enclone`, and will be described separately in more detail.

<img src="img/circles.svg" alt="circles" title="circles" />

___________________________________________________________________________________________________

### The `enclone` software

`enclone` is open-source, beta software.  Binary executables for Linux and Mac can be 
directly downloaded from this page, as can sample 10x Genomics datasets.  And then you're off and
running!  To use `enclone`, you need to know how to run command-line tools.  This is something that 
can be learned easily, particularly if you have a friend or colleague who can help you
get started.  You do not need to be able to program, or anything of that sort.

`enclone` is fast, typically responding in seconds.  It is intended as an experimental tool.
You can dynamically change your command line to select specific clonotypes and fields you wish
to see.

`enclone` is part of the [10x Genomics](https://www.10xgenomics.com/) immune 
profiling toolkit, including
[Cell Ranger and Loupe](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome), 
with which `enclone` will be integrated (later).

___________________________________________________________________________________________________

### How to understand `enclone` output

The example below shows how `enclone` prints out clonotypes.  This is something you'll need
to study in order to use `enclone` successfully.  `enclone` comes with extensive online 
documentation, and because you can easily play with the sample datasets, you can gradually
figure out how it all works.

<img src="img/enclone_annotated_example.svg" alt="enclone annotated example" title="enclone annotated example" /> By default, `enclone` prints clonotypes in this human-readable form.  You can also instruct
`enclone` to print clonotypes in machine-readable forms that are suitable for input to other
programs.

We are greatly interested in your feedback and ideas you may have to make `enclone` as useful
as possible.  You can write to us at enclone@10xgenomics.com.

___________________________________________________________________________________________________


### How to download and install `enclone`

1. Decide if you want to install `enclone` on a Linux or Mac computer.  Please let us know if
availability on other platforms is important to you.

2. Decide if you want the absolute latest version, or the most fully tested version.

3. Click one of these links to download the binary executables.

<b>TO DO: PUT LINKS HERE</b>

4. Move the file to an appropriate location.  For a desktop or laptop, an appropriate location 
might be `~/bin`.

<b>TO DO: link to pathing tutorial</b>

5. Change your path so that it includes this location.

6. Download the sample data.

The easiest way to do this is to type
```
git clone git@github.com:10XDev/enclone.git
```
and after doing that, `enclone/datasets` will be the directory containing the datasets
that are prepackaged with `enclone`.

7. Type `enclone help`, and read the terminal setup instructions there.

8. Type `enclone` to get to the main `enclone` help menu.

<b>Updates.</b> If you later choose to update `enclone`, you should update both the binary
executable and the datasets.  For the latter, you can simply type `git pull` inside the `enclone`
directory as created above.

___________________________________________________________________________________________________

### Questions

Please write to us at enclone@10xgenomics.com.  Please note that `enclone` is beta software
and that we may not be able to answer your question.

___________________________________________________________________________________________________

### How to compile the `enclone` software (for experts!)

Please see [compilation](COMPILE.md).
