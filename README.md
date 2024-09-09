# 2023-Transcriptome-hMSC

## About

Last version of repo before fully committing to DESeq2 pipeline.

## How to run

1.  Clone this repository

    -   ⚠️ cloning via RStudio highly recommended

2.  Add required data to specified folders under `/input`

    -   Mandatory:

        -   Data owned by Stem Cell and Neurogenesis Group - QUT

            -   `annotation`

            -   `cDNA`

            -   `cDNA_ncRNA_ENSEMBL`

        -   MSigDB gene sets

            -   `genesets`

3.  Open `2023-Transcriptome-hMSC.Rproj` in RStudio

    -   If this repo was cloned using RStudio, this project is already open.

4.  Source pipeline scripts at `/scripts/9_*` to perform analysis
