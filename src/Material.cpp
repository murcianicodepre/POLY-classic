#include "Material.h"

/*
    Material ~ Material class
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

Material::Material(float diff, float spec, float reflective, float refractive) : texture(NULL), bump(NULL), diff(diff), spec(spec), reflective(reflective), refractive(refractive) {}
void Material::loadTexture(const char* path){ texture = PolyRenderer::loadPNG(path); }
void Material::loadBump(const char* path){ bump = PolyRenderer::loadPNG(path); }