# Detecting illusory clonotype expansions

This page explains the origin of certain illusory clonotype expansions,
which might result in incorrect scientific conclusions, and shows how to detect them.

These illusory expansions are known to occur on occasion (see below for a demonstration), and we
hypothesize that they arise when an individual cell disintegrates or leaks, leaving fragments that 
seed multiple GEM partitions in the 10x system, thence yielding a clonotype which appears 
larger than its true size.

We believe that events of this type usually originate from plasma or plasmablast B cells.

Disintegration might occur during or after preparation of the sample.  One
way to document such an event would be to create two libraries from a single tube of cells.  If 
the clonotype is large and appears in only one of two libraries, one could be reasonably certain 
that a disintegration event occurred during or after cells were drawn from the tube.  This method 
could not be used to detect disintegration events occurring prior to that point.

Here we show that with the aid of gene expression data, illusory clonotype expansions can
generally be detected, even if only a single library was made.  The easier case would be a sample
consisting of pure B cells (for BCR).  The case where one has a mix of cell types (e.g. PBMCs) is 
more challenging because a GEM can contain both a B cell fragment (for BCR), plus a cell of a 
different type, and thus appear to have a normal level of gene expression, and no evidence of
mixing from the VDJ assay either.  We therefore focus on the case of samples that are mixed
cell types.

To that end, we show an example, using two libraries obtained from a single tube of PBMC cells, 
obtained from a healthy human donor.

```
enclone BCR=128037,128040 NCROSS
```

The `NCROSS` option instructs enclone to <i>not</i> filter out expanded clonotypes that appear
in only one dataset arising from the same sample (and which based on their sizes are highly
improbable).  Normally one would want this filtering, but these clonotypes are exactly what we
wish to see now!  Here is the top clonotype:

<img src="img/illusory1.png" alt="illusory1" title="illusory1" width=75% />

If we do not use the `NCROSS` option, but search for the clonotype using the heavy chain
CDR3 sequence, we see just one cell (the otherw having been filtered out)

```
enclone BCR=128037,128040 CDR3=CARGGTTTYFISW
```

<img src="img/illusory2.png" alt="illusory2" title="illusory2" width=75% />

Now suppose that both a VDJ and a GEX library have been made, as they have in this case.

```
enclone BCR=128040 GEX=127801 CDR3=CARGGTTTYFISW
```

<img src="img/illusory3.png" alt="illusory3" title="illusory3" width=75% />

Now we see less cells.  This is because the default behavior of enclone is to filter out
cells called by the VDJ pipeline that are not also called by the GEX pipeline.

Now we add the option `PER_CELL`, causing data for each cell to be displayed, and we also add two
fields to the display.  One is `gex`, the normalized count of gene expression UMIs,
and the other is a field `right`, that is more complicated.  We will also hide the onesie
(single chain) cells.

```
enclone BCR=128040 GEX=127801 CDR3=CARGGTTTYFISW PER_CELL LVARSP=gex,right MIN_CHAINS_EXACT=2
```

<img src="img/illusory4.png" alt="illusory4" title="illusory4" width=90% />

The field `right` is a measure of the extent to which cells having similar gene expression to a
 given putative B cell are 
themselves B cells (for BCR, or similarly for TCR).  In more detail, first let n be the number of 
VDJ cells that are also GEX cells.  Now for a given cell, we find the n GEX cells that are closest 
to it in PCA space, and report the percent of those that are also VDJ cells.  This is `right`.  The closer this number is to 100, the more the given cell looks like a typical B cell (or T cell, 
for TCR).  Conversely, a very low number makes the given cell appear suspect, although it is 
not <i>proof</i> of such.

The values of `right` vary considerably from dataset to dataset, requiring somewhat different
interpretation.  We show the distribution for this one dataset:

| right  | % of B cells  |
| -------| -------------:|
|  0-20  |  9.2          |
| 20-40  |  3.1          |
| 40-60  |  2.4          |
| 60-80  |  5.0          |
| 80-100 | 80.3          |

Thus the values of the cells in the reported clonotype are very low indeed, and almost all
highly suspect.  Probably the clonotype originated from a single cell, which broke up into one 
major piece (the one for barcode `CTGGTCTAGCTGCCCA-1`), and many smaller pieces.  These smaller 
pieces reside in GEMs that may or may not contain an actual intact cell.  In fact, 23 of the 
cells are detected as T cells (using TCR data).

We thus conclude in this case that the true clonotype probably consists of one cell.  Sometimes
one sees examples where there appear to be a few true cells, along with others that are not.
And sometimes one only sees only low UMI counts, likely corresponding to small fragments.
