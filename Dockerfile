FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL name="waldronlab/curatedMetagenomicAnalyses"

# For additional options in Jupyter-Server-Proxy
ENV RSESSION_PROXY_RSTUDIO_1_4="True"

# --- Optimized Caching ---
# Copy only the system installation script first. This layer will only be
# rebuilt if install.sh itself changes.
COPY scripts/install.sh /tmp/install.sh
RUN bash /tmp/install.sh \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# Copy the R installation script next. This layer is cached as long as
# install.R doesn't change.
COPY scripts/install.R /tmp/install.R
RUN R -f /tmp/install.R \
  && rm -f /tmp/install.R # Clean up the script after use

# By separating the installations, a change to install.R no longer
# requires re-running the system dependencies from install.sh.
# The original Dockerfile is also fine, this is just an optimization.

USER waldronlab
WORKDIR /home/waldronlab

EXPOSE 8888

# --- Add Healthcheck ---
# This command checks if the Jupyter Lab server is responding.
HEALTHCHECK CMD curl --fail http://localhost:8888/lab || exit 1

ENTRYPOINT ["/init"]

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--port=8888"]