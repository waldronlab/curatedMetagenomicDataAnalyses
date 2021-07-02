#!/usr/bin/with-contenv python3

from jupyter_server.auth import passwd
import os

password = os.getenv("PASSWORD", "waldronlab")

# Create a hash of the password and save it to the Jupyter Server Config
with open("/home/waldronlab/.jupyter/jupyter_server_config.py", "a") as f:
    f.write("c.ServerApp.password = '{}'".format(passwd(password)))
