cmake_minimum_required(VERSION 3.12)
project(PyBHJet)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify CMAKE_PREFIX_PATH for pybind11 
 #shoudl find the CMAKE_PREFIX_PATH for pybind11, can set it directly if it is not found. 

set(Python3_EXECUTABLE "pybhjet_venv/bin/python")

# set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};/path/to/pybind11/share/cmake/pybind11")
set(pybind11_DIR "pybhjet_venv/lib/python3.14/site-packages/pybind11/share/cmake/pybind11")

find_package(pybind11 REQUIRED)

find_package(GSL REQUIRED)

# environment variable for Python paths
if(NOT DEFINED ENV{PYTHONPATH})
    message(WARNING "PYTHONPATH environment variable is not set. Ensure pybind11 is accessible.")
endif()

# Add Kariba source files
set(KARIBA_SOURCES
    ../Kariba/BBody.cpp
    ../Kariba/Cyclosyn.cpp
    ../Kariba/Kappa.cpp
    ../Kariba/Mixed.cpp
    ../Kariba/Particles.cpp
    ../Kariba/ShSDisk.cpp
    ../Kariba/Compton.cpp
    ../Kariba/Bknpower.cpp
    ../Kariba/Powerlaw.cpp
    ../Kariba/Radiation.cpp
    ../Kariba/Thermal.cpp
)

set(PYBHJET_SOURCES
    python_interface.cpp
    bhjet_class.cpp
    jetmain.cpp
    utils.cpp
    jetpars.cpp
)

set(ALL_SOURCES ${PYBHJET_SOURCES} ${KARIBA_SOURCES})

# Add pybind11 module
pybind11_add_module(pybhjet ${ALL_SOURCES})

# Link libraries
target_link_libraries(pybhjet PRIVATE GSL::gsl GSL::gslcblas m pybind11::module)

# Include directories
target_include_directories(pybhjet PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/../Kariba
    $ENV{PYBIND11_INCLUDE_DIR}
    /opt/local/include
)

# Set library search paths
target_link_directories(pybhjet PRIVATE /opt/local/lib)