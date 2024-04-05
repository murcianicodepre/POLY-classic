#include "Camera.h"

/*
    Camera ~ world camera class
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

Camera::Camera(Vec3 ori, Vec3 rot, float fov) : ori(ori), rot(rot), fov(fov*ALPHA) {}
Ray Camera::rayTo(uint16_t x, uint16_t y){
    float aux = tanf(fov / 2.0f);
    float px = (2.0f * ((x + 0.5f)/WIDTH) - 1.0f) * aux * AR;
    float py = (1.0f - (2.0f * (y + 0.5f)/HEIGHT)) * aux;

    Vec3 dir = Vec3(px, py, 1.0f); dir.rotate(rot);
    return Ray(ori, dir.normalize());
}