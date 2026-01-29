#!/bin/bash

# Set the environment name
ENV_NAME="pybhjet_venv"

# Create the virtual environment
python3 -m venv $ENV_NAME

# Activate the virtual environment
source $ENV_NAME/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install required Python packages
pip install numpy pandas matplotlib jupyterlab pybind11 setuptools cmake gsl ipympl

# Check OS and install the correct compiler
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS - Use Clang
    brew install llvm@15
    export CC=$(brew --prefix llvm@15)/bin/clang
    export CXX=$(brew --prefix llvm@15)/bin/clang++
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux - Use GCC
    sudo apt update && sudo apt install -y gcc g++ make cmake
    export CC=gcc
    export CXX=g++
else
    echo "Unsupported OS. Please install a compatible compiler manually."
fi

# Install GSL
if [[ "$OSTYPE" == "darwin"* ]]; then
    brew install gsl
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    sudo apt install -y libgsl-dev
fi

# Print success message
echo "Virtual environment '$ENV_NAME' is ready."
echo "Activate it using: source $ENV_NAME/bin/activate"