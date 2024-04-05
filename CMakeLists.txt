# POLY classic ~ CMAKE file
# Diego Párraga Nicolás ~ diegojose.parragan@um.es

cmake_minimum_required(VERSION 3.16)
project(poly-classic LANGUAGES CXX)

# Search for packages
find_package(PNG REQUIRED)      # Read png textures and write png output
find_package(yaml-cpp)          # Parse scene script
find_package(OpenMP REQUIRED)   # OpenMP    
find_package(X11 REQUIRED)      # Poly GUI

set(CMAKE_INSTALL_PREFIX "/usr/local/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "..")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")    # Debug
set(CMAKE_CXX_STANDARD 17)

set(SOURCES                 # Set program sources
    src/main.cpp 
    src/Vec3.cpp 
    src/RGBA.cpp 
    src/Fragment.cpp 
    src/Ray.cpp 
    src/Camera.cpp
    src/Light.cpp
    src/Material.cpp
    src/Poly.cpp
    src/PolyRenderer.cpp 
)

add_executable(poly-classic ${SOURCES})

# Link libraries
target_link_libraries(poly-classic 
    PNG::PNG 
    ${YAML_CPP_LIBRARIES}
    OpenMP::OpenMP_CXX
    X11
)

# Add include folder
target_include_directories(poly-classic PRIVATE 
    "include"
    ${X11_INCLUDE_DIR}
)

install(TARGETS poly-classic    # Makefile output directory 
    RUNTIME DESTINATION "."
)