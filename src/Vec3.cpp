#include "Vec3.h"
#include "RGBA.h"

/*
    Vec3 ~ vec3 class for PolyRenderer
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

Vec3::Vec3() : x(0.0f), y(0.0f), z(0.0f) {}
Vec3::Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
Vec3::Vec3(float f) : x(f), y(f), z(f) {}
float Vec3::operator[] (const uint8_t i){ return i==0 ? x : (i==1 ? y : z); }
Vec3 Vec3::operator+ (Vec3 a) { return Vec3(x+a.x, y+a.y, z+a.z); }
Vec3 Vec3::operator- (Vec3 a) { return Vec3(x-a.x, y-a.y, z-a.z); }
Vec3 Vec3::operator* (Vec3 a) { return Vec3(x*a.x, y*a.y, z*a.z); }
Vec3 Vec3::operator/ (Vec3 a) { return Vec3(x/a.x, y/a.y, z/a.z); }
Vec3 Vec3::operator+ (float f) { return Vec3(x+f, y+f, z+f); }
Vec3 Vec3::operator- (float f) { return Vec3(x-f, y-f, z-f); }
Vec3 Vec3::operator* (float f) { return Vec3(x*f, y*f, z*f); }
Vec3 Vec3::operator/ (float f) { return Vec3(x/f, y/f, z/f); }
bool Vec3::operator== (float f) { return x==f && y==f & z==f; }
float Vec3::length() { return sqrtf(x*x + y*y + z*z); }
Vec3 Vec3::normalize() { float m = length(); return (m!=0.0f) ? Vec3(x/m, y/m, z/m) : Vec3(x,y,z); }
void Vec3::rotateX(float r){ 
    float a = r*ALPHA, y0 = y, z0 = z, cosa = cosf(a), sina = sinf(a);
    y = y0 * cosa - z0 * sina;
    z = y0 * sina + z0 * cosa;
}   
void Vec3::rotateY(float r){  
    float a = r*ALPHA, x0 = x, z0 = z, cosa = cosf(a), sina = sinf(a);
    x = x0 * cosa + z0 * sina;
    z = z0 * cosa - x0 * sina;
}  
void Vec3::rotateZ(float r){  
    float a = r*ALPHA, x0 = x, y0 = y, cosa = cosf(a), sina = sinf(a);
    x = x0 * cosa - y0 * sina;
    y = x0 * sina + y0 * cosa;
}  
void Vec3::rotate(Vec3 rot){
    this->rotateY(rot.y);
    this->rotateX(rot.x);
    this->rotateZ(rot.z);
}
float Vec3::dot(Vec3 a, Vec3 b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec3 Vec3::cross(Vec3 a, Vec3 b){ return Vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }
Vec3 Vec3::maxVec3(Vec3& a, Vec3& b){ return Vec3(a.x >= b.x ? a.x : b.x, a.y >= b.y ? a.y : b.y, a.z >= b.z ? a.z : b.z); } 
Vec3 Vec3::minVec3(Vec3& a, Vec3& b){ return Vec3(a.x <= b.x ? a.x : b.x, a.y <= b.y ? a.y : b.y, a.z <= b.z ? a.z : b.z); } 