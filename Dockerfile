FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL name="waldronlab/curatedMetagenomicAnalyses"

# For additional options in Jupyter-Server-Proxy
ENV RSESSION_PROXY_RSTUDIO_1_4="True"

# Copy both installation scripts
COPY scripts/install.sh /tmp/install.sh
COPY scripts/install.R /tmp/install.R

# Run both scripts in the same layer to ensure consistency
RUN bash /tmp/install.sh \
  && R -f /tmp/install.R \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/install.sh /tmp/install.R

USER waldronlab
WORKDIR /home/waldronlab

EXPOSE 8888

# --- Add Healthcheck ---
# This command checks if the Jupyter Lab server is responding.
HEALTHCHECK CMD curl --fail http://localhost:8888/lab || exit 1

ENTRYPOINT ["/init"]

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--port=8888"]