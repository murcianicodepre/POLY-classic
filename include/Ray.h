#ifndef __RAY_H__
#define __RAY_H__

/*
    Ray ~ poly-classic ray class 
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "POLY-classic.h"
#include "Vec3.h"

class Ray{
public:
    Vec3 ori, dir;
    uint32_t medium = 0u;   // By default, medium is "VOID_MAT"
    Ray();
    Ray(Vec3 ori, Vec3 dir);
    Vec3 point(float t);
};

#endif