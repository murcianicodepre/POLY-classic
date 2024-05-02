#ifndef __POLY_CLASSIC_H__
#define __POLY_CLASSIC_H__

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
#include <vector>
#include <string>
#include <cstring>
#include <csignal>
#include <optional>
#include "png.h"
#include "yaml-cpp/yaml.h"
#include "yaml-cpp/exceptions.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

typedef uint32_t uint;

// Render pipeline defs
constexpr uint16_t WIDTH = 1280;
constexpr uint16_t HEIGHT = 960;
constexpr uint16_t TEXTURE_SIZE = 1024;
constexpr uint8_t TILE_SIZE = 8;
constexpr float AR = 1.33333f;
constexpr uint8_t MAX_RAY_BOUNCES = 255u;
constexpr bool BLENDER = true;


// Individual tri rendering flags
constexpr uint8_t DISABLE_RENDERING = 0x01u;
constexpr uint8_t DISABLE_TEXTURES = 0x02u;
constexpr uint8_t DISABLE_SHADING = 0x04u;
constexpr uint8_t DISABLE_BUMP = 0x08u;
constexpr uint8_t DISABLE_TRANSPARENCY = 0x10u;
constexpr uint8_t DISABLE_SHADOWS = 0x20u;
constexpr uint8_t DISABLE_REFLECTIONS = 0x40u;
constexpr uint8_t DISABLE_REFRACTIONS = 0x80u;

// Rendering pipeline flags
constexpr uint8_t DISABLE_FAST_INTERSECTION_SHADER = 0x01u;
constexpr uint8_t FLAT_SHADING = 0x02u;

// Math defs
#define PI 3.14159265f
#define ALPHA 0.017453292f
#define EPSILON 1e-6f

// Forward class declaration
class Vec3;
class RGBA;
class Fragment;
class Ray;
class Tri;
class Vertex;
class Material;
struct Hit;

// Shared classes
static const char* getCpu();
static void printIntro();
static void polyMsg(std::string msg);
Vec3 parseVec3(YAML::Node node);
RGBA parseColor(YAML::Node node);
uint8_t parseFlag(YAML::Node node);

#endif