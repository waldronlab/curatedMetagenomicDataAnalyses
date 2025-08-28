#!/usr/bin/bash

REPOSITORY="curatedMetagenomicAnalyses"
PASSWORD=${PASSWORD:=waldronlab}

# Install JupyterHub dependencies
apt-get update
apt-get install -y --no-install-recommends npm nodejs libzmq3-dev
npm install -g configurable-http-proxy
rm -rf ~/.npm

RUN python3 -m venv /env
ENV PATH="/env/bin:$PATH"
RUN pip install --upgrade pip

python3 -m pip install jupyterlab jupyter-server-proxy jupyter-rsession-proxy jupyter_client jupyter_core
  
# Install Python analysis dependencies
apt-get install -y --no-install-recommends python3-skbio python-skbio-doc
python3 -m pip install pandas \
  setuptools \
  numpy \
  python-dateutil \
  pytz \
  statsmodels \
  matplotlib \
  seaborn \
  num2words

# Add default user
useradd -m -G staff,rstudio waldronlab 
echo "waldronlab:$PASSWORD" | chpasswd

# Allow user to write to RStudio DB
chmod -R g+w /var/lib/rstudio-server
chown -R :rstudio /var/lib/rstudio-server

# Remove s6 Rstudio service
rm -rf /etc/services.d/rstudio

# Create Jupyter configuration file
cd /home/waldronlab
su -c "jupyter server --generate-config" waldronlab

# Remove password and token requirement
cat <<EOT >> /home/waldronlab/.jupyter/jupyter_server_config.py
c.NotebookApp.password = ""
c.NotebookApp.token = ""
c.NotebookApp.ip = "*"
EOT

# Clone repository and change permissions waldronlab
git clone --depth 1 https://github.com/waldronlab/$REPOSITORY.git /home/waldronlab/$REPOSITORY
chown -R waldronlab:waldronlab /home/waldronlab/$REPOSITORY

# Copy README.md and adjust links for analyses paths
cat /home/waldronlab/$REPOSITORY/README.md | sed "s/(vignettes/($REPOSITORY\/vignettes/g" >> /home/waldronlab/README.md
