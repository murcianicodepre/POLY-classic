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
#define AR 1.333333333f
#define MAX_REFLECTIONS 8u
#define MAX_REFRACTIONS 8u
#define BLENDER true
constexpr uint8_t DISABLE_RENDERING = 0x01u;
constexpr uint8_t DISABLE_TEXTURES = 0x02u;
constexpr uint8_t DISABLE_SHADING = 0x04u;
constexpr uint8_t DISABLE_BUMP = 0x08u;

// Math defs
#define PI 3.14159265f
#define ALPHA 0.017453292f
#define EPSILON 1e-6f

// Incomplete classes 
class V3f;
class Fragment;
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
    PolyRenderer();
    ~PolyRenderer();
    RGBA* frame;                    // Rendering buffer
    Camera* cam;                    // Scene camera
    std::vector<Tri> tris;          // Tri data
    std::vector<Material> mats;     // Material data
    std::vector<Light> lights;      // Scene light data
    bool loadScene(const char* scene);
    bool render(uint threads);
    void save(const char* path);
    RGBA compute_pixel(uint x, uint y);
    bool intersection_shader(Ray&, Hit&, uint32_t DISCARD = 0u);
    Fragment fragment_shader(Hit&);
    Fragment texture_shader(Hit& hit);
    V3f bump_shader(Hit& hit);
    Fragment reflection_shader(Ray&, Hit&, uint);
    Fragment refraction_shader(Ray&, Hit&, uint);
};

#endif