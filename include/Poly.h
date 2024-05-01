#ifndef __POLY_H__
#define __POLY_H__

/*
    Poly ~ Vertex, Hit, Tri and Poly class headers
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "POLY-classic.h"
#include "Vec3.h"
#include "Ray.h"
#include "PolyRenderer.h"

/* 
    Vertex class, contains:
    - World coordinates
    - Normal vector
    - UV texture coordinates
*/
class Vertex{
public:
    Vec3 xyz, normal;
    float u, v;
    Vertex() : xyz(), normal(), u(), v() {}
    Vertex(Vec3 xyz, Vec3 normal, float u, float v);
    void move(Vec3 m);
    void scale(Vec3 s);
    void scale(float s);
    void rotate(Vec3 r);
    void rotateX(float r);
    void rotateY(float r);
    void rotateZ(float r);
};

/*
    Hit record struct, stores:
    - Tri global identifier
    - Hit point, geometric and phong normals
    - Texture coordinates of the hit
*/
struct Hit{
    uint32_t tri = 0u;
    Vec3 point, normal, phong;
    float u = 0.0f, v = 0.0f;
    Ray ray;
    bool valid = false;
};

// Tri class
class Tri {
public:
    Vertex a, b, c;
    uint16_t mat, poly;
    Vec3 centroid;
    uint8_t flags;
    Tri(Vertex a, Vertex b, Vertex c, uint16_t poly, uint16_t mat, uint8_t flags);
    bool intersect(Ray ray, Hit& hit);
    void move(Vec3 m);
    void scale(float s);
    void scale(Vec3 s);
    void rotate(Vec3 r);
    void rotateX(float r);
    void rotateY(float r);
    void rotateZ(float r);
    float min(uint8_t axis);
    float max(uint8_t axis);
};

// Poly class
class Poly {
public:
    std::vector<Tri> tris;
    Poly(const char* path, uint16_t polyId, uint16_t mat, uint8_t flags);
    void move(Vec3 m);
    void scale(Vec3 s);
    void scale(float s);
    void rotate(Vec3 r);
    void rotateX(float r);
    void rotateY(float r);
    void rotateZ(float r);
};

#endif