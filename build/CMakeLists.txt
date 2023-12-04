# POLY classic ~ CMAKE file
# Diego Párraga Nicolás ~ diegojose.parragan@um.es

cmake_minimum_required(VERSION 3.16)
project(poly-classic LANGUAGES CXX)

# Search for packages
find_package(PNG REQUIRED)  # Read png textures and write png output
find_package(yaml-cpp)      # Parse scene script

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "..")
set(CMAKE_CXX_STANDARD 14)

set(SOURCES                 # Set program sources
    ../src/main.cpp 
    ../src/poly-classic.cpp 
    ../include/poly-classic.h
)

# Link libraries
target_link_libraries(poly-classic 
    PNG::PNG 
    ${YAML_CPP_LIBRARIES}
)

# Add include folder
target_include_directories(poly-classic PRIVATE 
    "../include"
)

install(TARGETS poly-classic    # Makefile output directory 
    RUNTIME DESTINATION "."
)