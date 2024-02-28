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
#include <tuple>
#include <string>
#include <cstring>
#include "png.h"
#include "yaml-cpp/yaml.h"
#include "yaml-cpp/exceptions.h"

// Render pipeline defs
#define WIDTH 1280
#define HEIGHT 960
#define TEXTURE_SIZE 1024
#define FOV 60.0f
#define AR 1.333333333f
#define BLENDER true
#define DISABLE_RENDERING 0x01u
#define DISABLE_SHADING 0x02u
#define DISABLE_TEXTURES 0x04u
#define DISABLE_BUMP 0x8u

// Math defs
#define PI 3.14159265f
#define ALPHA 0.017453292f
#define EPSILON 1e-4f

// Incomplete classes 
class RGBA;
class Camera;
class Light;
class Tri;
class Material;
class Vertex;
class Hit;
class Ray;

typedef unsigned int uint;


// Other functions
const char* getCpuCode();
void printIntro();
void polymsg(std::string msg);

// Main renderer class
class PolyRenderer{
public:
    RGBA* frame;                    // Rendering buffer
    Camera* cam;                    // Scene camera
    std::vector<Tri> tris;          // Tri data
    std::vector<Material> mats;     // Material data
    std::vector<Light> lights;      // Scene light data
    PolyRenderer();
    ~PolyRenderer();
    bool loadScene(const char* scene);
    bool render(uint threads);
    void save(const char* path);
    RGBA compute_pixel(uint x, uint y);
    bool intersection_shader(Ray&, Hit&);
    RGBA fragment_shader(Hit&);
    RGBA texture_shader(Hit& hit);
};

#endif