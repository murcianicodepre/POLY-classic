#include "Ray.h"

/*
    Ray ~ poly-classic ray header 
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

Ray::Ray() : ori(), dir() {}
Ray::Ray(Vec3 ori, Vec3 dir) : ori(ori), dir(dir) {}
Vec3 Ray::point(float t){ return dir*t + ori; }
float Ray::getT(Vec3 point){ return ((point-ori) / dir).x; }