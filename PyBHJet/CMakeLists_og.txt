cmake_minimum_required(VERSION 3.12)
project(PyBHJet)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify CMAKE_PREFIX_PATH for pybind11 
 #shoudl find the CMAKE_PREFIX_PATH for pybind11, can set it directly if it is not found. 

set(Python3_EXECUTABLE "${CMAKE_SOURCE_DIR}/pybhjet_venv/bin/python")

## forcing it to be found, but only looking in the venv (so not great)
execute_process( 
  COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
  OUTPUT_VARIABLE pybind11_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(pybind11_DIR "${pybind11_dir}")
find_package(pybind11 CONFIG REQUIRED)
find_package(GSL REQUIRED)

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
target_link_libraries(pybhjet PRIVATE GSL::gsl GSL::gslcblas m)

# Include directories
target_include_directories(pybhjet PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/../Kariba
    /opt/local/include
)

# Set library search paths
target_link_directories(pybhjet PRIVATE /opt/local/lib)