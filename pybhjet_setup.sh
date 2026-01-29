#!/bin/bash
echo "Building PyBHJet..."
mkdir -p build && cd build
cmake ..
make
echo "Build complete!"