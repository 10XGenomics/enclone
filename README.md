<a name="readme" style="display:block; position:relative; top:-150px;"></a>
# enclone

<img align="left" src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/cf/Construction_workers_not_wearing_fall_protection_equipment.jpg/320px-Construction_workers_not_wearing_fall_protection_equipment.jpg" alt="dangerous constuction zone" title="dangerous construction zone" />

<br>
<br>
<br>

Watch your step!  This site is under construction and undergoing initial testing.  Please do not 
enter unless you've been directed here.

(Photo credit: NIOSH/Wikipedia)

<br>
<br>
<br>
<br>

<b>enclone</b> is a computational tool for studying the adaptive immune system in humans and 
other vertebrate species.  It takes as input data from 
<b>[10x Genomics](https://www.10xgenomics.com/)</b>, providing captured RNA 
sequences for B and T cell receptors in single cells.  It organizes these cells into groups 
(clonotypes) arising from the same progenitor and compactly displays each clonotype along 
with its salient features, including mutated amino acids.

enclone and this page are designed for immunologists, but anyone can download and experiment 
with it.

<b>Background:</b> when you get sick, your body mounts an immune response by selectively amplifying 
immune cells and mutations within these selected cells.  enclone allows you to see the history of 
single immune cells within a biological sample (such as a blood draw or biopsy).  This history 
reflects how the cognate receptors of these cells evolved in response to antigens, including 
viruses, bacteria, and tumors.
___________________________________________________________________________________________________

## The mission of enclone

If you have a sample and have generated single-cell VDJ data from the B or T cells within it, you
have the power to fully understand the nature of the receptors for those cells, because you have
in hand the sequences for each of their chains.  

You should be able to directly see this biology, without aid from a computational expert.  To 
that end, we start with this simple goal:

<img align="left" src="img/mission.svg" alt="mission" title="mission" />
The <i>finding</i> part is 
algorithmically challenging: it is very easy to mistakenly put unrelated cells in the same 
clonotype, or to "pollute" a clonotype with chains that do not belong in it.
The <i>displaying</i> part is also challenging: we are unaware of other tools that do so.

<br>
<br>
<br>

To make sure we are using terminology in the same way, the following diagram shows what a 
_clonotype_ is for B cells.  The same applies for T cells, but things are simpler because T cells 
do not have somatic hypermutation.
<img src="img/what_is_a_clonotype.svg" alt="what is a clonotype" title="what is a clonotype" />

Each cell in a clonotype is typically represented by two or three chains.  Such information can
only be obtained from _single cell_ data!  From such data, clonotypes can be computationally
approximated, with high accuracy (see below).  The method we use for this is described briefly in 
the online documentation for enclone, and will be described separately in more detail.

<img src="img/performance.svg" alt="performance" title="performance" />

___________________________________________________________________________________________________

<a name="software" style="display:block; position:relative; top:-150px;"></a>
## The enclone software

enclone is open-source, beta software.  Binary executables for Linux and Mac can be 
directly downloaded from this page, as can sample 10x Genomics datasets.  And then you're off and
running!  To use enclone, you need to know how to run command-line tools.  This is something that 
can be learned easily, particularly if you have a friend or colleague who can help you
get started.  You do not need to be able to program, or anything of that sort.

enclone is fast, typically responding in seconds (if run on a single dataset).  It is intended 
as an exploratory tool.  You can dynamically change your command line to select specific 
clonotypes and fields you wish to see.  You can run enclone on a laptop or a desktop or a server.

enclone is part of the [10x Genomics](https://www.10xgenomics.com/) immune 
profiling toolkit, including
[Cell Ranger and Loupe](https://support.10xgenomics.com/single-cell-vdj), 
with which enclone will be integrated (later).
___________________________________________________________________________________________________

<a name="download" style="display:block; position:relative; top:-150px;"></a>
## Installing enclone

<b>(these instructions are not yet functional)</b>

<b>1.  Open terminal window.</b>  Open a terminal window on your Linux or Mac computer. Please let us
know if availability on other platforms is important to you.

<b>2.  Download enclone.</b>  Type the following to download the enclone executable:
```
mkdir -p ~/bin; cd ~/bin
wget https://github.com/10XGenomics/enclone/releases/download/latest/linux/enclone
or on a Mac
wget https://github.com/10XGenomics/enclone/releases/download/latest/mac/enclone
```
This gets you the absolute latest version of enclone.  You can repeat this step if you ever
want to update.  At a later date, there will also be separately numbered releases that have passed 
a more extensive set of tests.

It is not necessary to compile enclone, unless you want to contribute
to the enclone codebase.  Please see [compilation](COMPILE.md).

<b>3.  Download test data.</b>  Type the following to download the enclone test datasets 
(plus source code, because it's easier to fetch everything):
```
cd
svn export https://github.com/10XGenomics/enclone/trunk enclone
```
(See [here](fetching_test_datasets.md#readme) if this doesn't work for you.)  At this point 
`~/enclone/datasets` will contain the datasets that are prepackaged with enclone.  If you 
subsequently want to update this, delete the directory and repeat the command.

<b>4.  Update your path.</b>  Edit your shell initialization file to add `:~/bin` to `PATH`.  Ask a colleague for help
if needed.  Close and reopen your terminal window to refresh your path.  Then you're good to go!
___________________________________________________________________________________________________

## Running enclone

Running enclone can be as simple as typing e.g. 
```
enclone BCR=/home/my_name/experiment_123
```
where the path is where your Cell Ranger outputs live, but there are many options to learn
about.  For example, if you want to combine many datasets, you can do that, but you probably
need to provide a metadata file that describes the datasets.  You can find most of the enclone
documentation within its online menus.  To get started you should:

1. Type `enclone help`, to make sure your terminal window works for `enclone`.

2. Type `enclone` to get to the main enclone help menu.

The concatenated help pages are also
[here](https://htmlpreview.github.io/?https://github.com/10XGenomics/enclone/blob/master/src/help.all.html) 
<b>[BROKEN LINK]</b>.  We may expand this out in the future to show the separate pages.
___________________________________________________________________________________________________

## Understanding enclone output

The example below shows how enclone prints out clonotypes.  This is something you'll need
to study in order to use enclone successfully.  enclone comes with extensive online 
documentation, and because you can easily play with the sample datasets, you can gradually
figure out how it all works.

<img src="img/enclone_annotated_example.svg" alt="enclone annotated example" title="enclone annotated example" /> 

Notice the compression in two directions.  Vertically, rather than showing one line for every cell,
we group cells into a single line if they have identical VDJ transcripts.  Horizontally, rather than
showing all transcript positions, we only show "interesting" positions.  This is a flexible 
concept, and what we show by default are all positions exhibiting a difference from the reference
and all positions in the CDR3.

The same exact output would be obtained by typing
```
enclone PRE=~/enclone/datasets BCR=123085 CDR3=CQQRSNWPPSITF
```
provided that you put test data in the location indicated under download instructions.  Otherwise
you would need to change the value of `PRE`.  The directory `123085` is in the directory
`~/enclone/datasets` and contains some files from a Cell Ranger run, obtained from a human 
ovarian cancer sample.

The argument `CDR3=CQQRSNWPPSITF` causes enclone to display only clonotypes in which the given
CDR3 sequence occurs.  Many other filters are provided.  In the absence of filters, all clonotypes
are shown.  Clonotypes are shown from largest to smallest, and the output is automatically paged,
so you can scroll through it.

By default, enclone prints clonotypes in this human-readable form.  You can also instruct enclone 
to print clonotypes in machine-readable forms that are suitable for input to other programs.
___________________________________________________________________________________________________

## Combining multiomic data

Gene expression and feature barcode data can be displayed simultaneously with VDJ data.  For
example, here we add columms for the same clonotype, showing the median number of UMIs detected
for, all genes, a particular gene, and a particular antibody:

<img src="img/clonotype_with_gex.png" alt="clonotype with gex" title="clonotype with gex" /> 

To obtain this, we added the extra arguments
```GEX=123749 LVARSP=gex,IGHV3-49_g,CD19_ab```
to the previous command.  The `GEX` part points to the directory containing gene expression and
feature barcode data.  The `LVARSP` part defines the additional columns to be displayed.

Other types of data can be brought in via featuring barcoding.  For example, response to 
multiple antigens can be measured using [LIBRA-seq](https://www.ncbi.nlm.nih.gov/pubmed/31787378)
and these data displayed as additional columns.
___________________________________________________________________________________________________

<a name="honeycomb" style="display:block; position:relative; top:-150px;"></a>
## Visualizing many clonotypes at once

<img align="left" src="img/clono.svg" alt="honeycomb plot" title="honeycomb plot" />

You can select clonotypes in enclone and then display them using a "honeycomb" plot.

In this instance, pre- and post-vaccination samples were collected from four individuals,
many datasets were generated for each sample, and these were combined in a single call
to enclone.  Clonotypes containing at least ten cells are shown.  The plot was generated by adding
```
MIN_CELLS=10 PLOT="clono.svg,pre->blue,post->red 
LEGEND=blue,"pre-vaccination cell",
       red,"post-vaccination cell"
```
to the enclone command line, yielding the image shown here as the file `clono.svg`.

<br><br><br><br>
___________________________________________________________________________________________________

## So what is enclone good for?

There are many ways to use 10x data to study immune biology.  Thus in the previous section, the 
red clonotypes may represent responses to antigens in the vaccine.  But for a given 
disease <i>e.g.</i> COVID-19, one may not yet have a vaccine, and indeed the vaccine may be
the goal!  A different approach is needed.

One such approach is to <b> identify patient and survivor B cell
clonotypes that expand in response to infectious disease</b>. These define antibodies that can be 
used to design passive or active vaccines.  Additional power is added by mapping 
antigen specificity to multiple antigens directly via feature barcoding
([LIBRA-seq](https://www.ncbi.nlm.nih.gov/pubmed/31787378)).

These data are easy to display in enclone!  You can directly select candidates for 
vaccine or therapeutic development by picking large clonotypes with high antigen counts and single 
or multiple antigen specifities.

We are actively working on further functionality that will make this process even more effective.
___________________________________________________________________________________________________

## Questions

Please ask us questions!  We are greatly interested in your feedback and ideas you may have to 
make enclone as useful as possible.  You can write to us at enclone@10xgenomics.com.
Please note that enclone is beta software and thus a work in progress.  We are actively
making many changes and may be unable to respond promptly to your particular request.
