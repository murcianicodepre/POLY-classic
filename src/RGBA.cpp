#include "RGBA.h"
#include "Fragment.h"
#include "Vec3.h"

using namespace std;

/*
    RGBA ~ rgba class for PolyRenderer
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

RGBA::RGBA() : r(0u), g(0u), b(0u), a(0u) {}
RGBA::RGBA(uint8_t r, uint8_t g, uint8_t b, uint8_t a) : r(r), g(g), b(b), a(a) {}
RGBA::RGBA(Fragment frag){
    Fragment aux = frag;
    aux.clamp(1.0f, 0.0f);
    
    r = static_cast<uint8_t>(frag.r*255.0f);
    g = static_cast<uint8_t>(frag.g*255.0f);
    b = static_cast<uint8_t>(frag.b*255.0f);
    a = static_cast<uint8_t>(frag.a*255.0f);
}
Vec3 RGBA::toVec3(){
    return Vec3(static_cast<float>(r)/255.0f, static_cast<float>(g)/255.0f, static_cast<float>(b)/255.0f);
}