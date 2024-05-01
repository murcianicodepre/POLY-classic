#ifndef __VEC3_H__
#define __VEC3_H__

/*
    Vec3 ~ vec3 header for PolyRenderer
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "POLY-classic.h"

using namespace std;

class Vec3{
public:
    float x,y,z;
    Vec3();
    Vec3(float x, float y, float z);
    Vec3(float f);
    float operator[] (const uint8_t i);
    Vec3 operator+ (Vec3 a);
    Vec3 operator- (Vec3 a);
    Vec3 operator* (Vec3 a);
    Vec3 operator+ (float f);
    Vec3 operator- (float f);
    Vec3 operator* (float f);
    bool operator== (float f);
    float length();
    Vec3 normalize();
    void rotateX(float r);
    void rotateY(float r);
    void rotateZ(float r);
    void rotate(Vec3 rot);

    static float dot(Vec3 a, Vec3 b);
    static Vec3 cross(Vec3 a, Vec3 b);
    static Vec3 maxVec3(Vec3& a, Vec3& b);
    static Vec3 minVec3(Vec3& a, Vec3& b);
};

#endif 