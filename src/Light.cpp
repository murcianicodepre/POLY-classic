#include "Light.h"

/*
    Light ~ light class
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

Ambient::Ambient(RGBA color, float intensity) : Light(color, intensity) {}
LightType Ambient::type() { return LightType::Ambient; }

PointLight::PointLight(Vec3 pos, RGBA color, float intensity) : Light(color, intensity), pos(pos) {}
LightType PointLight::type(){ return LightType::Point; }

DirectionalLight::DirectionalLight(Vec3 dir, RGBA color, float intensity) : Light(color, intensity), dir(dir) {}
LightType DirectionalLight::type(){ return LightType::Direction; }