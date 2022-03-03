#/usr/bin/env bash
# ------------------------------
# Script name: autopilot_master.sh
#
# Script purpose: run Change-O and SCOPer on 10x input single cell data
# Author: Dr. Wyatt J. McDonnell
# Acknowledgments: Dr. David B. Jaffe
# Date created: 2022-01-26
# Copyright (c) 10x Genomics, Inc. 2022
# Email: wyatt.mcdonnell@10xgenomics.com
# ------------------------------
#
# Notes:
# This script makes use of the Immcantation Docker image. The latest version of
# the Docker container at the time of writing was 4.3.0. This script does the
# following:
#
# 1) Generates mandatory inputs for Immcantation using enclone's
#    build_immcantation_inputs binary. The binary for build_immcantation_inputs
#    can be built as part of building enclone. This generates the following
#    files:
#    .
#    ├── filtered_contig_annotations.csv
#    └── filtered_contig.fasta
#
# 2) Reports information about how many barcodes go into Immcantation
# 3) Re-annotates the contigs using Change-O and IMGT references found in the
#    Immcantation Docker container. This generates the following files:
#    .
#    ├── filtered_contig_db-pass.tsv
#    ├── filtered_contig_heavy_productive-F.tsv
#    ├── filtered_contig_heavy_productive-T.tsv
#    ├── filtered_contig_igblast.fmt7
#    ├── filtered_contig_light_productive-F.tsv
#    ├── filtered_contig_light_productive-T.tsv
#    ├── logs
#    |   ├── pipeline-10x.err
#    |   └── pipeline-10x.log
#    └── temp_files.tar.gz
#
# 4) Runs SCOPer's filtering routines and its hierarchicalClones() function
#    to group VDJ sequences into clonotypes, and writes the clonotypes out to a
#    TSV file. This is achieved by using Rscript to run the autopilot.R script
#    also provided in this repo. This generates the following files:
#    .
#    ├── scoper_clones.tsv
#    ├── scoper_filter.tsv
#    └── scoper.log
#
# 5) When the "pipeline" finishes, the final directory structure should look
#    like this:
#    .
#    ├── autopilot.R
#    ├── autopilot_test.sh
#    ├── filtered_contig_annotations.csv
#    ├── filtered_contig_db-pass.tsv
#    ├── filtered_contig.fasta
#    ├── filtered_contig_heavy_productive-F.tsv
#    ├── filtered_contig_heavy_productive-T.tsv
#    ├── filtered_contig_igblast.fmt7
#    ├── filtered_contig_light_productive-F.tsv
#    ├── filtered_contig_light_productive-T.tsv
#    ├── logs
#    │   ├── pipeline-10x.err
#    │   └── pipeline-10x.log
#    ├── scoper_clones.tsv
#    ├── scoper_filter.tsv
#    ├── scoper.log
#    └── temp_files.tar.gz
#
# This script takes as input the following files:
# 1. filtered_contig_heavy_productive-T.tsv
# 2. filtered_contig_light_productive-T.tsv
# 3. the name of a file to output productive contigs and associated data
# 4. the name of a file to output clonotyped data
#
# R is finicky and can break. We apologize if this happens! The version of
# Rscript used in the Docker container is as follows for 4.3.0 of the
# Immcantation Docker container:
# `R scripting front-end version 4.0.5 (2021-03-31)`
#
# Rscript in the container can be found at /usr/bin/Rscript.
# ------------------------------start=`date +%s`
echo "Starting workflow from $PWD..."
echo "Building Immcantation inputs..."
# N.B. you can remove ./ if build_immcantation_inputs is on your $PATH
./build_immcantation_inputs @test
echo "Immcantation inputs generated!"
# Print the number of unique cell barcodes in the contig annotations
nbarcodesin=$(awk -F ',' '{print $1}' filtered_contig_annotations.csv | sort | uniq | wc -l)
echo "There are $nbarcodesin barcodes in the contig annotations."
date
echo "Running Change-O..."
# N.B. you may wish to use fewer than 32 threads on your machine
docker run -v $PWD:/data:z immcantation/suite:4.3.0 changeo-10x -s data/filtered_contig.fasta -a data/filtered_contig_annotations.csv -o data -g human -t ig -p 32 -o data
date
echo "Change-O complete!"
date
echo "Running SCOPer..."
docker run -v $PWD:/data:z immcantation/suite:4.3.0 Rscript data/autopilot.R data/filtered_contig_heavy_productive-T.tsv data/filtered_contig_light_productive-T.tsv scoper_filter.tsv scoper_clones.tsv
# Report information about # of barcodes at each step of the analysis
nbarcodesfilter=$(awk -F '\t' '{print $49}' scoper_filter.tsv | sort | uniq | wc -l)
nbarcodesout=$(awk -F '\t' '{print $49}' scoper_clones.tsv | sort | uniq | wc -l)
echo "$nbarcodesfilter out of $nbarcodesin barcodes were retained prior to clonotyping."
echo "$nbarcodesout out of $nbarcodesin barcodes were retained after Change-O and SCOPer."
end=`date +%s`
runtime=$((end-start))
echo "Runtime was $runtime seconds."
