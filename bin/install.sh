#! /usr/bin/bash

USER="waldronlab"
REPOSITORY="curatedMetagenomicAnalyses"

# Install JupyterHub dependencies

apt-get update \
  && apt-get install -y --no-install-recommends npm nodejs libzmq3-dev \
  && npm install -g configurable-http-proxy \
  && rm -rf ~/.npm \
  && python3 -m pip install jupyterhub notebook jupyterlab
  # jupyter-server-proxy

# Install Python analyses dependencies

python3 -m pip install pandas \
  setuptools \
  numpy \
  python-dateutil \
  pytz \
  statsmodels

# Add default user

useradd -m -G staff $USER \
  && echo '$USER:$USER' | chpasswd

# Clone repository as waldronlab

su -c $(git clone https://github.com/waldronlab/$REPOSITORY.git \
  /home/$USER/$REPOSITORY) $USER
