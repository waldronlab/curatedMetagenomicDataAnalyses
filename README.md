| Resource               | link                                                                                                         |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------- |
| Data                   | |[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17567900.svg)](https://doi.org/10.5281/zenodo.17567900)      |
| Analyses               | [This Repo](https://github.com/waldronlab/curatedMetagenomicDataAnalyses/tree/main/cMD3_paper_analyses)           |
| CuratedMetagenomicData | [View in Bioconductor](https://bioconductor.org/packages/release/data/experiment/html/curatedMetagenomicData.html)|

# curatedMetagenomicDataAnalyses

This repository provides biologically relevant analyses using the [curatedMetagenomicData](https://bioconductor.org/packages/curatedMetagenomicData/) package, both using R/Bioconductor and using Python. You can run both R and Python analyses locally in the provided Docker container, or on the Cloud for free.

## Running locally using Docker

### Requirements

You need [Docker](https://docs.docker.com/get-docker/).

### Getting Started

First build the image:

    docker pull ghcr.io/waldronlab/curatedmetagenomicanalyses:latest

Then run a container based on the image with your password:

    docker run -d -p 80:8888 --name cma \
      waldronlab/curatedmetagenomicanalyses

Visit `localhost/lab` in your browser.

## Running locally without Docker

Start with an installation of the current version of Bioconductor (see https://bioconductor.org/install/). Older versions probably will not work. 
Installation directly from GitHub requires first installing the `remotes` package, then:
```r
BiocManager::install("waldronlab/curatedMetagenomicDataAnalyses", dependencies = TRUE)
```

## Analyses

### R Vignettes

* [Create datasets for machine learning](https://waldronlab.io/curatedMetagenomicDataAnalyses/articles/MLdatasets.html)
* [Exploration of the liver cirrhosis dataset](https://waldronlab.io/curatedMetagenomicDataAnalyses/articles/explorecirrhosis.html)
* [Select all colorectal cancer patients and create a cohort table, calculate prevalence of all species found in their stool microbiomes and create a dynamic searchable html table](https://waldronlab.io/curatedMetagenomicDataAnalyses/articles/identify_CRC_species.html)
* [Meta-analysis of age-related microbial species using cMD3](https://waldronlab.io/curatedMetagenomicDataAnalyses/articles/Age_metaanalysis_vignette.html)
* [Meta-analysis of sex-related microbial species using cMD3](https://waldronlab.io/curatedMetagenomicDataAnalyses/articles/Sex_metaanalysis_vignette.html)
* [NUI Galway Metagenomics Workshop](https://waldronlab.io/curatedMetagenomicDataAnalyses/articles/NUI-Galway-Metagenomics-Workshop.html)

### Python Notebooks

* [Sex-related differences in the human microbiome using cMD3 and Python3](https://github.com/waldronlab/curatedMetagenomicDataAnalyses/blob/main/vignettes/sexContrastMicrobiomeAnalysis.ipynb)

### Supplementary Materials 

* [Installing Python dependencies in Linux](https://github.com/waldronlab/curatedMetagenomicDataAnalyses/blob/main/vignettes/installation.ipynb) (Python notebook)
