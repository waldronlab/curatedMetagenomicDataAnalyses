FROM bioconductor/bioconductor_docker:RELEASE_3_13

LABEL name="waldronlab/curatedMetagenomicAnalyses"

COPY bin /tmp

# Install Jupyterhub
RUN bash /tmp/install.sh

# Install IRkernel for R notebooks and curatedMetagenomicAnalyses dependencies
RUN R -f /tmp/install.R

# WORKDIR /home/waldronlab

CMD ["jupyterhub"]
