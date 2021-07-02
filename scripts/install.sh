#!/usr/bin/bash

REPOSITORY="curatedMetagenomicAnalyses"
PASSWORD=${PASSWORD:=waldronlab}

# Install JupyterHub dependencies
apt-get update
apt-get install -y --no-install-recommends npm nodejs libzmq3-dev
npm install -g configurable-http-proxy
rm -rf ~/.npm
python3 -m pip install jupyterlab

# Install repositories to facilitate Jupyterlab proxying RStudio
for repo in jupyter-server-proxy jupyter-rsession-proxy
do
  git clone https://github.com/jupyterhub/$repo.git /tmp/$repo
  python3 -m pip install /tmp/$repo
done
rm -rf /tmp/jupyter*

# Install Python analysis dependencies
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

# Create Jupyter configuration file
cd /home/waldronlab
su -c "jupyter server --generate-config" waldronlab

# Move script to s6 initialization script path
mv /tmp/userhashconf.py /etc/cont-init.d/
chmod 755 /etc/cont-init.d/userhashconf.py

# Clone repository and change permissions waldronlab
git clone https://github.com/waldronlab/$REPOSITORY.git /home/waldronlab/$REPOSITORY
chown -R waldronlab:waldronlab /home/waldronlab/$REPOSITORY
