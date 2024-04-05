#ifndef __RGBA_H__
#define __RGBA_H__

/*
    RGBA ~ rgba header for PolyRenderer
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "POLY-classic.h"

class RGBA{
public:
    uint8_t r,g,b,a;
    RGBA();
    RGBA(uint8_t r, uint8_t g, uint8_t b, uint8_t a);
    RGBA(Fragment frag);
    Vec3 toVec3();
};

#endif