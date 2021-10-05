library(BiocManager)

# Install R kernel for R Jupyter
BiocManager::install("rzmq")
BiocManager::install("repr")
BiocManager::install("IRkernel")
BiocManager::install("IRdisplay")
IRkernel::installspec(user = FALSE);

# Install dependencies for curatedMetagenomicAnalyses
BiocManager::install("waldronlab/curatedMetagenomicAnalyses", dependencies = TRUE)
BiocManager::install("waldronlab/curatedMetagenomicDataCLI", dependencies = TRUE)
curatedMetagenomicDataCLI::install()
