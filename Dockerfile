# POLY-classic dockerfile ~ diegojose.parragan@um.es
FROM alpine:latest

# Update distro and install dependencies
RUN apk update && apk add --no-cache gcc g++ cmake make libpng-dev openmp git

# First step, clone yaml-cpp and compile the static library
WORKDIR /
RUN git clone https://github.com/jbeder/yaml-cpp.git
WORKDIR /yaml-cpp 
RUN mkdir build && cd build && cmake .. && make && make install

# Second step, copy source inside /app and compile poly-classic
WORKDIR /app
COPY CMakeLists.txt CMakeLists.txt
COPY include include
COPY src src
RUN mkdir build
WORKDIR /app/build
RUN cmake .. && make && make install

# Last remove source code and set working dir
WORKDIR /poly-classic
RUN rm -fr /app

# Image is ready to run poly-classic passing the input .poly scene description
# Example: docker run -v PATH-TO-SCENE-DATA:/poly-classic poly-classic SCENE.POLY [NTHREADS]