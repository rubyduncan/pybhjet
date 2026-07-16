#!/bin/bash

# Set the environment name
ENV_NAME="pybhjet_venv"

# Create the virtual environment (make sure that on helios you load something higher than 3.6)
python3 -m venv $ENV_NAME

# Activate the virtual environment
source $ENV_NAME/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install required Python packages
pip install numpy pandas matplotlib astropy jupyterlab pybind11 setuptools cmake gsl ipympl ipykernel

#this version is working at present, should keep to these for now (feb. 2026)
pip install "astromodels==2.4.3"
pip install "threeml==2.4.3"

#can be used to ignore the XSPEC verison because it will try (and fail) to find it 
# env -u HEADAS -u ASTRO_XSPEC_VERSION -u XSPEC_INC_PATH python -m pip install --no-cache-dir astromodels

# Print success message
echo "Virtual environment '$ENV_NAME' is ready."
echo "Activate it using: source $ENV_NAME/bin/activate"