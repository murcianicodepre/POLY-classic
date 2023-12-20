#ifndef _POLY_CLASSIC_
#define _POLY_CLASSIC_

/*
    POLY classic ~ main header
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <omp.h>
#include <cmath>
#include <list>
#include <vector>
#include <string>
#include <cstring>
#include "png.h"
#include "yaml-cpp/yaml.h"
#include "yaml-cpp/exceptions.h"

// Render pipeline defs
#define WIDTH 1024
#define HEIGHT 768
#define TEXTURE_SIZE 1024
#define FOV 60.0f
#define AR 1.333333333f

// Math defs
#define PI 3.14159265f
#define ALPHA 0.017453292f
#define EPSILON 1e-4f

// Incomplete classes 
class RGBA;
class Camera;
class Tri;

// Main renderer class
class PolyRenderer{
public:
    RGBA* frame;
    Camera* cam;

    PolyRenderer();
    ~PolyRenderer();
    bool loadScene(const char* path);
    bool render(uint threads);
    void save(const char* path);
};

#endif