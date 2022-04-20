#!/usr/bin/Rscript --vanilla
# ------------------------------
# Script name: autopilot.R
#
# Script purpose: run SCOPer on 10x input single cell data
# Author: Dr. Wyatt J. McDonnell
# Acknowledgments: Drs. Kenneth Hoehn and Hailong Meng
# Date created: 2022-01-26
# Copyright (c) 10x Genomics, Inc. 2022
# Email: wyatt.mcdonnell@10xgenomics.com
# ------------------------------
#
# Notes:
# This script makes use of the R language found in the Immcantation Docker
# image. The latest version of the Docker container at the time of writing was
# 4.3.0. This script processes cell barcode-associated contigs that have passed
# the changeo-10x pipeline and uses SCOPer's hierarchical clustering method to
# group the contigs (and cells) into clonotypes. We have observed difficulties
# with successful parallelization of this pipeline using the authors' `nproc`
# argument and so this variable is set to nproc=1. It may behave differently on
# your machine. We were unable to run the spectral clustering mode of SCOPer on
# several datasets >1M single cells even allowing 7 days for runtime; you can
# use the spectralClones() function in place of hierarchicalClones() if you wish
# to try this yourself.
#
# This script requires the following files to be present in the directory:
# 1. filtered_contig_heavy_productive-T.tsv
# 2. filtered_contig_light_productive-T.tsv
# 3. the name of a file to output productive contigs and associated data
# 4. the name of a file to output clonotyped data
#
# Note that these files are automatically generated as part of
# autopilot_master.sh, also provided in this repository, and which executes
# this script.
#
# R is finicky and can break. We apologize if this happens! The version of
# Rscript used in the Docker container is as follows for 4.3.0 of the
# Immcantation Docker container:
# `R scripting front-end version 4.0.5 (2021-03-31)`
#
# Rscript within the Immcantation Docker container can be found at
# /usr/bin/Rscript.
# ------------------------------

### load libraries
library(alakazam)
library(shazam)
library(dplyr)
library(tidyr)
library(data.table)
library(scoper)

###
setwd('/data')

### memory monitoring
memStart <- gc(reset=TRUE)

### initialize tidy dataframe
data = tibble()

### read in data
#### reads in Change-O "productive" heavy chain sequences
h = readChangeoDb('filtered_contig_heavy_productive-T.tsv')
#### reads in Change-O "productive" light chain sequences
l = readChangeoDb('filtered_contig_light_productive-T.tsv')
#### joins heavy and light rows; requires column symmetry
comb = bind_rows(h,l)
data = comb

### filter to 'productive' contigs; used in case wrong contig passed in
data = filter(data, productive == TRUE)
### writes out sequences prior to grouping
writeChangeoDb(data, 'scoper_filter.tsv')

### find clones
multi_heavy = table(filter(data, locus=="IGH")$cell_id)
multi_heavy_cells = names(multi_heavy)[multi_heavy > 1]
### remove multi-heavy-chain cells
data = filter(data, !cell_id %in% multi_heavy_cells)
### note that we used thresholds of [0.05, 0.10, 0.15, 0.20] in the paper
clones = hierarchicalClones(data, threshold=0.1, cell_id='cell_id', only_heavy=FALSE, split_light=TRUE, verbose=TRUE, log='scoper.log', summarize_clones=FALSE)
### write out clonotyped sequences
writeChangeoDb(clones, 'scoper_clones.tsv')
gc()
memEnd <- gc()
message(paste('SCOPer complete!'))
message(paste(memEnd[12] - memStart[12], "megabytes were allocated during the SCOPer part of the workflow."))
