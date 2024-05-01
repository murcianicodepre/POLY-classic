#ifndef __FRAGMENT_H__
#define __FRAGMENT_H_

/*
    Fragment ~ fragment header for PolyRenderer
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "POLY-classic.h"

class Fragment{
public:
    float r,g,b,a;
    Fragment();
    Fragment(float r, float g, float b, float a);
    Fragment(Vec3 rgb, float a);
    Fragment(RGBA color);
    void clamp(float max, float min);
    Fragment operator+ (Fragment frag);
    Fragment operator- (Fragment frag);
    Fragment operator* (Fragment frag);
    Fragment operator* (Vec3 vec);
    Fragment operator+ (float f);
    Fragment operator- (float f);
    Fragment operator* (float f);
};

#endif