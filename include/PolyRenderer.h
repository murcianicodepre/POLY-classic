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
    RGBA* frame;
    Camera* cam;
    std::vector<Tri> tris;
    std::vector<Material> mats;
    std::vector<std::unique_ptr<Light>> lights;
    BVHNode* bvh;

    // Main program functions
    bool loadScene(const char* scene);
    bool render(uint8_t threads, bool ENABLE_RENDERING_WINDOW);
    void save(const char* path);

    // Acceleration struct
    uint32_t nextNode = 1, * triIdx;
    void buildBVH();
    void updateNodeBounds(uint32_t nodeId);
    void subdivide(uint32_t nodeId);
    void intersectBVH(Ray& ray, Hit& hit, uint32_t nodeId, uint16_t discard = 0x0000u);

    // Rendering pipeline shaders
    RGBA compute_pixel(uint16_t x, uint16_t y);
    bool intersection_shader(Ray& ray, Hit& hit, uint16_t discard = 0x0000u);
    Fragment blinn_phong_shading(Hit& hit);
    Fragment flat_shading(Hit& hit);
    //Fragment fragment_shader(Hit& hit);
    Fragment fragment_shader(Hit&);
    Fragment texture_mapping(Hit& hit);
    Vec3 bump_mapping(Hit& hit);
    Fragment reflection_shader(Hit&, uint8_t);
    Fragment refraction_shader(Hit&, uint8_t);
    Fragment raytracing_shader(Hit&, uint8_t, uint8_t);

    // Other renderer functions
    static RGBA* loadPNG(const char* path);
    static void savePNG(const char* path, RGBA* texture);
    static const char* getCpu();
    static void printIntro();
    static void polyMsg(std::string msg);
    static Vec3 parseVec3(YAML::Node node);
    static RGBA parseColor(YAML::Node node);
    static uint16_t parseFlags(YAML::Node node);

    // Debug flags
    uint16_t debug;
}; 



#endif 