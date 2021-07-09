FROM bioconductor/bioconductor_docker:RELEASE_3_13

LABEL name="waldronlab/curatedMetagenomicAnalyses"

# For additional options in Jupyter-Server-Proxy
ENV RSESSION_PROXY_RSTUDIO_1_4="True"

COPY scripts /tmp

# Install Jupyterhub, IRkernel, and curatedMetagenomicAnalyses dependencies
RUN bash /tmp/install.sh \
  && R -f /tmp/install.R \
  && rm -rf /tmp/install* \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

USER waldronlab

WORKDIR /home/waldronlab

EXPOSE 8888

ENTRYPOINT ["/init"]

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--port=8888"]
