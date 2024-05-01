#ifndef __CAMERA_H__
#define __CAMERA_H__

/*
    Camera ~ world camera class header
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "POLY-classic.h"
#include "Vec3.h"
#include "Ray.h"

class Camera{
public:
    Vec3 ori, lookAt;
    float fov;
    Camera(Vec3 ori, Vec3 lookAt, float fov);
    Ray rayTo(uint16_t x, uint16_t y);
};

#endif