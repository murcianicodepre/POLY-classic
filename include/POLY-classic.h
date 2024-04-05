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
#include <thread>
#include <cmath>
#include <list>
#include <vector>
#include <tuple>
#include <string>
#include <cstring>
#include <optional>
#include <bitset>
#include <memory>
#include "png.h"
#include "yaml-cpp/yaml.h"
#include "yaml-cpp/exceptions.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>

// Render pipeline defs
#define WIDTH 1280
#define HEIGHT 960
#define TEXTURE_SIZE 1024
#define TILE_SIZE 8
#define AR 1.333333333f
#define MAX_REFLECTIONS 8u
#define MAX_REFRACTIONS 8u
#define BLENDER true

// Individual tri rendering flags
constexpr uint8_t DISABLE_RENDERING = 0x01u;
constexpr uint8_t DISABLE_TEXTURES = 0x02u;
constexpr uint8_t DISABLE_SHADING = 0x04u;
constexpr uint8_t DISABLE_BUMP = 0x08u;
constexpr uint8_t DISABLE_TRANSPARENCY = 0x10u;
constexpr uint8_t DISABLE_SHADOWS = 0x20u;

// Debug flags
constexpr uint8_t DRAW_AABB = 0x01u;
constexpr uint8_t RENDER_TRIS = 0x02u;
constexpr uint8_t RENDER_NORMALS = 0x04u;

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