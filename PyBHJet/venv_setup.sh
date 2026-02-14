#!/bin/bash

ENV_NAME="pybhjet_venv"

python3 -m venv $ENV_NAME
source $ENV_NAME/bin/activate
pip install --upgrade pip

# Install required Python packages - just for pybhjet 
pip install numpy pandas matplotlib astropy jupyterlab pybind11 setuptools cmake gsl ipympl 

# if threeML will be used; uncomment this line as well, and make sure your astromodels install is complete first
#pip install astromodels 
#pip install threeml

# Can be used to ignore the XSPEC verison because it may try (and fail) to find it if using threeml
# env -u HEADAS -u ASTRO_XSPEC_VERSION -u XSPEC_INC_PATH python -m pip install --no-cache-dir astromodels

# Print success message
echo "Virtual environment '$ENV_NAME' is ready."
echo "Activate it using: source $ENV_NAME/bin/activate"
