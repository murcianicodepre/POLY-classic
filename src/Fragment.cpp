#include "Fragment.h"
#include "RGBA.h"
#include "Vec3.h"

/*
    Fragment ~ fragment class for PolyRenderer
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

void Fragment::clamp(float max, float min){
    r = r > max ? max : (r<min ? min : r);
    g = g > max ? max : (g<min ? min : g);
    b = b > max ? max : (b<min ? min : b);
    a = a > max ? max : (a<min ? min : a);
}
Fragment::Fragment() : r(0.0f), g(0.0f), b(0.0f), a(0.0f) {}
Fragment::Fragment(float r, float g, float b, float a) : r(r), g(g), b(b), a(a) {}
Fragment::Fragment(Vec3 rgb, float a) : r(rgb.x), g(rgb.y), b(rgb.z), a(a) {}
Fragment::Fragment(RGBA color) : 
    r(static_cast<float>(color.r)/255.0f),
    g(static_cast<float>(color.g)/255.0f),
    b(static_cast<float>(color.b)/255.0f),
    a(static_cast<float>(color.a)/255.0f)
{}
Fragment Fragment::operator+ (Fragment frag){ Fragment out = Fragment(r+frag.r, g+frag.g, b+frag.b, a+frag.a); out.clamp(1.0f, 0.0f); return out; }
Fragment Fragment::operator- (Fragment frag){ Fragment out = Fragment(r-frag.r, g-frag.g, b-frag.b, a-frag.a); out.clamp(1.0f, 0.0f); return out; }
Fragment Fragment::operator* (Fragment frag){ Fragment out = Fragment(r*frag.r, g*frag.g, b*frag.b, a*frag.a); out.clamp(1.0f, 0.0f); return out; }
Fragment Fragment::operator+ (float f){ Fragment out = Fragment(r+f, g+f, b+f, a+f); out.clamp(1.0f, 0.0f); return out; }
Fragment Fragment::operator- (float f){ Fragment out = Fragment(r-f, g-f, b-f, a-f); out.clamp(1.0f, 0.0f); return out; }
Fragment Fragment::operator* (float f){ Fragment out = Fragment(r*f, g*f, b*f, a*f); out.clamp(1.0f, 0.0f); return out; }