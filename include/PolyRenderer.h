#ifndef __POLYRENDERER_H__
#define __POLYRENDERER_H__

/*
    PolyRenderer ~ main poly-classic renderer header 
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "POLY-classic.h"
#include "Vec3.h"
#include "Ray.h"
#include "RGBA.h"
#include "Fragment.h"
#include "Light.h"
#include "Camera.h"
#include "Poly.h"

using namespace std;

class BVHNode;

class PolyRenderer{
public:
    PolyRenderer();
    ~PolyRenderer();
    RGBA* _frame;
    Camera* _cam;
    std::vector<Tri> _tris;
    std::vector<Material> _mats;
    std::vector<std::unique_ptr<Light>> _lights;
    BVHNode* _bvh = nullptr;

    // Main program functions
    bool loadScene(const char* scene);
    bool render(uint8_t threads);
    void save(const char* path);

    // Acceleration struct
    uint32_t _nextNode = 1, * _triIdx = nullptr, _splits = 128u;
    void buildBVH();
    void updateNodeBounds(uint32_t nodeId);
    void subdivide(uint32_t nodeId);
    void intersectBVH(Ray& ray, Hit& hit, uint32_t nodeId, uint16_t discard = 0x0000u);
    float EvaluateSAH(BVHNode&, uint8_t, float);
    float getBestSplit(BVHNode&, uint8_t&, float&);
    float getNodeCost(BVHNode&);

    // Rendering pipeline functions
    RGBA compute_pixel(uint16_t, uint16_t);
    bool intersection_shader(Ray&, Hit&, uint8_t discard = 0x00u);
    Fragment blinn_phong_shading(Hit&, uint8_t flags = 0x00u);
    Fragment flat_shading(Hit&);
    Fragment compute_fragment(Hit&, uint8_t flags = 0x00u);
    Fragment texture_mapping(Hit&);
    Vec3 bump_mapping(Hit&);
    Fragment fragment_shader(Hit&, uint8_t);
    Fragment compute_reflection(Hit&, uint8_t);
    Fragment compute_refraction(Hit&, uint8_t);

    // Other renderer functions
    static RGBA* loadPNG(const char* path);
    static void savePNG(const char* path, RGBA* texture);
    static const char* getCpu();
    static void printIntro();
    static void polyMsg(std::string msg);
    static Vec3 parseVec3(YAML::Node node);
    static RGBA parseColor(YAML::Node node);
    static uint16_t parseFlags(YAML::Node node, bool print = false);

    // Global flags
    uint16_t _global = 0x00000000u;
}; 

#endif 