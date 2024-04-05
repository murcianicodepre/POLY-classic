#ifndef __LIGHT_H__
#define __LIGHT_H__

/*
    Light ~ light class header
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "POLY-classic.h"
#include "Vec3.h"
#include "RGBA.h"

enum class LightType { Ambient, Point, Direction };

// Base light class
class Light{
public:
    RGBA color;
    float intensity;
    Light(RGBA color, float intensity) : color(color), intensity(intensity) {}
    virtual LightType type() = 0;
};

// Ambient light source
class Ambient : public Light{
public:
    Ambient(RGBA color, float intensity);
    LightType type() override;
};

// Point light source
class PointLight : public Light{
public:
    Vec3 pos;
    PointLight(Vec3 pos, RGBA color, float intensity);
    LightType type() override;
};

// Directional light source
class DirectionalLight : public Light{
public:
    Vec3 dir;
    DirectionalLight(Vec3 dir, RGBA color, float intensity);
    LightType type() override;
};

#endif 