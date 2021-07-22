library(BiocManager)

# Install Rkernel for R Jupyter
install.packages(c("rzmq", "repr", "IRkernel", "IRdisplay"))
IRkernel::installspec(user = FALSE);

# Install dependencies for curatedMetagenomicAnalyses
BiocManager::install("waldronlab/curatedMetagenomicAnalyses",
                     dependencies = TRUE)
BiocManager::install("waldronlab/curatedMetagenomicDataCLI",
                     dependencies = TRUE)
curatedMetagenomicDataCLI::install()
