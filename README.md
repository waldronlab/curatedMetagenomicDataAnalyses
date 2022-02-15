# Curated Metagenomic Analyses

This repository provides biologically relevant analyses using the [curatedMetagenomicData](https://bioconductor.org/packages/curatedMetagenomicData/) package, both using R/Bioconductor and using Python. You can run both R and Python analyses locally in the provided Docker container, or on the Cloud for free.

## Running on the Cloud for free

A machine with all dependencies, code from this repository, and Jupyterlab (with R and Python3) and RStudio running is available at http://app.orchestra.cancerdatasci.org/ (search for the Curated Metagenomic Analyses workshop). You can use these machines for up to 8 hours at a time.

## Running locally using Docker

### Requirements

You need [Docker](https://docs.docker.com/get-docker/).

### Getting Started

First build the image:

    docker build -t "waldronlab/curatedmetagenomicanalyses" .

Then run a container based on the image with your password:

    docker run -d -p 80:8888 --name cma \
      waldronlab/curatedmetagenomicanalyses

Visit `localhost` in your browser.

## Analyses

### R Vignettes

* [Create datasets for machine learning](https://waldronlab.io/curatedMetagenomicAnalyses/articles/MLdatasets.html)
* [Exploration of the liver cirrhosis dataset](https://waldronlab.io/curatedMetagenomicAnalyses/articles/explorecirrhosis.html)
* [Meta-analysis of age-related microbial species using cMD3](https://waldronlab.io/curatedMetagenomicAnalyses/articles/Age_metaanalysis_vignette.html)
* [Meta-analysis of sex-related microbial species using cMD3](https://waldronlab.io/curatedMetagenomicAnalyses/articles/Sex_metaanalysis_vignette.html)
* [NUI Galway Metagenomics Workshop](https://github.com/waldronlab/curatedMetagenomicAnalyses/blob/vignettes/NUI-Galway-Metagenomics-Workshop.Rmd)

### Python Notebooks

* [Sex-related differences in the human microbiome using cMD3 and Python3](https://github.com/waldronlab/curatedMetagenomicAnalyses/blob/main/vignettes/sexContrastMicrobiomeAnalysis.ipynb)

### Supplementary Materials 

* [Installing Python dependencies in Linux](https://github.com/waldronlab/curatedMetagenomicAnalyses/blob/main/vignettes/installation.ipynb) (Python notebook)
