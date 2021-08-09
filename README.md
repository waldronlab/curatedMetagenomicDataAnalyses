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

* [Create datasets for machine learning](vignettes/MLdatasets.Rmd)
* [Exploration of the liver cirrhosis dataset](vignettes/explorecirrhosis.Rmd)
* [Meta-analysis of age-related microbial species using cMD3](vignettes/Age_metaanalysis_vignette.Rmd)
* [Meta-analysis of sex-related microbial species using cMD3](vignettes/Sex_metaanalysis_vignette.Rmd)

### Python Notebooks

* [Sex-related differences in the human microbiome using cMD3 and Python3](vignettes/sexContrastMicrobiomeAnalysis.ipynb)

### Supplementary Materials 

* [Installing Python dependencies in Linux](vignettes/installation.ipynb) (Python notebook)
