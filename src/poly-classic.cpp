/*
    POLY-classic ~ main header class
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "../include/poly-classic.h"
using namespace std;

std::vector<std::string> msgvector;

/* Cpu arch code */
const char* getCpuCode(){
    char vendor[13];
    uint eax = 0, ebx, ecx, edx;

    __asm__ __volatile__(
        "cpuid"
        : "=b"(ebx), "=c"(ecx), "=d"(edx)
        : "a"(eax)
    );

    memcpy(vendor, &ebx, 4);
    memcpy(vendor + 4, &edx, 4);
    memcpy(vendor + 8, &ecx, 4);
    vendor[12] = '\0';
    return vendor[0] == 'G' ? "\e[1;94mx86-64" : "\e[1;91mamd64";
}

/* Print intro */
void printIntro(){
    system("clear");
    printf(" \e[1;91m▄▄▄   \e[92m▄▄   \e[94m▄  \e[95m▄   ▄ \n");
    printf(" \e[91m█  █ \e[92m█  █  \e[94m█   \e[95m█ █  \e[93m ▄▄ ▄   ▄   ▄▄  ▄▄ ▄  ▄▄\n");
    printf(" \e[91m█▀▀  \e[92m█  █  \e[94m█    \e[95m█  \e[93m █   █  █▄█ ▀▄  ▀▄  █ █\n");
    printf(" \e[91m█    \e[92m▀▄▄▀  \e[94m█▄▄  \e[95m█   \e[93m▀▄▄ █▄ █ █ ▄▄▀ ▄▄▀ █ ▀▄▄\n");
    printf("\e[91m         - diegojose.parragan@um.es -\n\e[0m\n");
}

/* Print parser msg and save it to buffer for reprinting */
void polymsg(string msg){
    printf("%s", msg.c_str());
    msgvector.push_back(msg);
}

/* Main vector class */
class V3f {
public:
    float x, y, z;
    V3f() : x(0.0f), y(0.0f), z(0.0f) {}
    V3f(float x, float y, float z) : x(x), y(y), z(z) {}
    V3f operator+ (V3f a) { return V3f(x+a.x, y+a.y, z+a.z); }
    V3f operator- (V3f a) { return V3f(x-a.x, y-a.y, z-a.z); }
    V3f operator* (V3f a) { return V3f(x*a.x, y*a.y, z*a.z); }
    V3f operator+ (float f) { return V3f(x+f, y+f, z+f); }
    V3f operator- (float f) { return V3f(x-f, y-f, z-f); }
    V3f operator* (float f) { return V3f(x*f, y*f, z*f); }
    float length() { return sqrtf(x*x + y*y + z*z); }
    V3f normalize() { float m = length(); return V3f(x/m, y/m, z/m); }
    void rotateX(float r){ 
        float a = r*ALPHA, y0 = y, z0 = z, cosa = cosf(a), sina = sinf(a);
        y = y0 * cosa - z0 * sina;
        z = y0 * sina + z0 * cosa;
    }   
    void rotateY(float r){  
        float a = r*ALPHA, x0 = x, z0 = z, cosa = cosf(a), sina = sinf(a);
        x = x0 * cosa + z0 * sina;
        z = z0 * cosa - x0 * sina;
    }  
    void rotateZ(float r){  
        float a = r*ALPHA, x0 = x, y0 = y, cosa = cosf(a), sina = sinf(a);
        x = x0 * cosa - y0 * sina;
        y = x0 * sina + y0 * cosa;
    }  
    string toString(){ return "{" + to_string(x) + ", " + to_string(y) + ", " + to_string(z) + "}"; }
};
float dot(V3f a, V3f b){ return a.x*b.x + a.y*b.y + a.z*b.z; }  
V3f cross(V3f a, V3f b){ return V3f(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }  

/* RGBA pixel class */
class RGBA{
public:
    uint16_t r, g, b, a;
    RGBA() : r(0u), g(0u), b(0u), a(0u) {}
    RGBA(uint16_t r, uint16_t g, uint16_t b, uint16_t a) : r(r), g(g), b(b), a(a) {}
    V3f asV3f(){ return V3f(static_cast<float>(r)/255.0f, static_cast<float>(g)/255.0f, static_cast<float>(b)/255.0f); }
    RGBA(V3f v) : r(static_cast<uint16_t>(v.x*255.0f)), g(static_cast<uint16_t>(v.y*255.0f)), b(static_cast<uint16_t>(v.z*255.0f)), a(255) {}
    void clamp255(){ r = r>255 ? 255 : r; g = g>255 ? 255 : g; b = b>255 ? 255 : b; a = a>255 ? 255 : a; }
    RGBA operator+ (RGBA color) { RGBA out = RGBA(r+color.r, g+color.g, b+color.b, a+color.a); out.clamp255(); return out; }
    RGBA operator- (RGBA color) { RGBA out = RGBA(r-color.r, g-color.g, b-color.b, a-color.a); out.clamp255(); return out; }
    RGBA operator* (float f) { RGBA out = RGBA(static_cast<uint16_t>(r*f), static_cast<uint16_t>(g*f), static_cast<uint16_t>(b*f), a); out.clamp255(); return out; }
    RGBA operator* (V3f v) { RGBA out = RGBA(static_cast<uint16_t>(r*v.x), static_cast<uint16_t>(g*v.y), static_cast<uint16_t>(b*v.z), a); out.clamp255(); return out; }
};

/* Fragment class. RGBA for graphics, range of values: 0.0f ~ 1.0f */
class Fragment{
public:
    float r,g,b,a;
    Fragment() : r(0.0f), g(0.0f), b(0.0f), a(0.0f) {}
    Fragment(float r, float g, float b, float a) : r(r), g(g), b(b), a(a) {}
    Fragment(V3f vector, float a) : r(vector.x), g(vector.y), b(vector.z), a(a) {}
    Fragment(RGBA color) : r(static_cast<float>(color.r)/255.0f), g(static_cast<float>(color.g)/255.0f), b(static_cast<float>(color.b)/255.0f), a(static_cast<float>(color.a)/255.0f) {}
    void clamp255(){ 
        r = r>1.0f ? 1.0f : (r<0.0f ? 0.0f : r); g = g>1.0f ? 1.0f : (g<0.0f ? 0.0f : g); 
        b = b>1.0f ? 1.0f : (b<0.0f ? 0.0f : b); a = a>1.0f ? 1.0f : (a<0.0f ? 0.0f : a); 
    }
    Fragment operator+ (Fragment frag) { Fragment out = Fragment(r+frag.r, g+frag.g, b+frag.b, a+frag.a); out.clamp255(); return out; }
    Fragment operator- (Fragment frag) { Fragment out = Fragment(r-frag.r, g-frag.g, b-frag.b, a-frag.a); out.clamp255(); return out; }
    Fragment operator* (Fragment frag) { Fragment out = Fragment(r*frag.r, g*frag.g, b*frag.b, a*frag.a); out.clamp255(); return out; }
    Fragment operator* (float f) { Fragment out = Fragment(r*f, g*f, b*f, a*f); out.clamp255(); return out; }
    RGBA toRGBA(){ return RGBA(static_cast<uint8_t>(r*255.0f), static_cast<uint8_t>(g*255.0f), static_cast<uint8_t>(b*255.0f), static_cast<uint8_t>(a*255.0f)); }
};

/* Load png as RGBA texture */
RGBA* loadPNG(const char* path){
    FILE* input = fopen(path, "rb");
        if(!input){ polymsg("\e[1;91m  err loading '" + string(path) + "': file could not be opened!\n\e[0m"); exit(EXIT_FAILURE); }

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if(!png){ 
            fclose(input);
            polymsg("\e[1;91m  err texture '" + string(path) + "' read struct could not be loaded!\n\e[0m"); exit(EXIT_FAILURE);
        }

    png_infop info = png_create_info_struct(png);
        if(!info){
            fclose(input);
            png_destroy_read_struct(&png, NULL, NULL);
            polymsg("\e[1;91m  err texture '" + string(path) + "' info struct could not be loaded!\e[0m\n"); exit(EXIT_FAILURE);
        }

    png_init_io(png, input);
    png_read_info(png, info);

    int w = png_get_image_width(png, info), h = png_get_image_height(png, info);
    png_byte color_type = png_get_color_type(png, info);
    png_byte bit_depth = png_get_bit_depth(png, info);

    if(bit_depth==16) png_set_strip_16(png);
    if(color_type==PNG_COLOR_TYPE_RGB) png_set_palette_to_rgb(png);
    if(color_type==PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png);
    if(png_get_valid(png, info, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);
    if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_PALETTE) 
        png_set_filler(png, 0xff, PNG_FILLER_AFTER);
    if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) png_set_gray_to_rgb(png);

    png_read_update_info(png, info);

    RGBA* texture = (RGBA*) malloc(sizeof(RGBA) * TEXTURE_SIZE * TEXTURE_SIZE);
    vector<png_bytep> row_p(h);
    for(int y=0; y<h; y++)
        row_p[y] = new png_byte[png_get_rowbytes(png, info)];

    png_read_image(png, row_p.data());

    for(int y=0; y<h; y++)
        for(int x=0; x<w; x++){
            png_bytep pixel = &(row_p[y][x * 4]); 
            texture[y * w + x].r = pixel[0];
            texture[y * w + x].g = pixel[1];
            texture[y * w + x].b = pixel[2];
            texture[y * w + x].a = pixel[3];
        }

    for(int y=0; y<h; y++)
        delete[] row_p[y];

    png_destroy_read_struct(&png, &info, NULL);

    if(fclose(input)==-1){ polymsg("\e[1;91m  err closing '" + string(path) + "'\e[0m\n"); exit(EXIT_FAILURE); }
    polymsg("\e[1;93m  loaded texture '" + string(path) + "'\e[0m\n");
    
    return texture;
}

/* Save RGBA* to png */
void savePNG(const char* path, RGBA* texture){
    FILE* output = fopen(path, "wb");
        if(!output){ polymsg("\e[1;91m  err saving to '" + string(path) + "': file could not be opened!\n\e[0m"); exit(EXIT_FAILURE); }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if(!png){ 
            fclose(output);
            polymsg("\e[1;91m  err saving render\n\e[0m");
            exit(EXIT_FAILURE);
        }

    png_infop info = png_create_info_struct(png);
        if(!info){ 
            fclose(output);
            png_destroy_write_struct(&png, NULL);
            polymsg("\e[1;91m  err saving render\e[0m\n");
            exit(EXIT_FAILURE);
        }
    
    png_init_io(png, output);
    png_set_IHDR(png, info, WIDTH, HEIGHT, 8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);

    vector<png_byte> row(WIDTH * 4);
    for(int y=0; y<HEIGHT; y++){
        for(int x=0; x<WIDTH; x++){
            row[x * 4] = static_cast<png_byte>(texture[y * WIDTH + x].r);               // r
            row[x * 4 + 1] = static_cast<png_byte>(texture[y * WIDTH + x].g);          // g
            row[x * 4 + 2] = static_cast<png_byte>(texture[y * WIDTH + x].b);         // b    
            row[x * 4 + 3] = static_cast<png_byte>(texture[y * WIDTH + x].a);        // a
        }
        png_write_row(png, row.data());
    }
    png_write_end(png, NULL);
    png_destroy_write_struct(&png, &info);

    if(fclose(output)==-1){ polymsg("\e[1;91m  err closing '" + string(path) + "'\e[0m\n"); exit(EXIT_FAILURE); }
}

/* Ray class */
class Ray {
public:
    V3f ori, dir;
    Ray() : ori(), dir() {}
    Ray(V3f ori, V3f dir) : ori(ori), dir(dir) {}
    V3f point(float t){ return dir*t + ori; }
};

/* Camera class */
class Camera {
public:
    V3f ori;
    float fov, rx, ry, rz;
    Camera(V3f ori, float fov) : ori(ori), fov(fov*ALPHA), rx(0.0f), ry(0.0f), rz(0.0f) {}
    Ray rayTo(int x, int y){
        float aux = tanf(fov / 2.0f);
        float px = (2.0f * ((x + 0.5f)/WIDTH) - 1.0f) * aux * AR;
        float py = (1.0f - (2.0f * (y + 0.5f)/HEIGHT)) * aux;

        V3f dir = V3f(px, py, 1.0f);
            dir.rotateZ(rz);
            dir.rotateY(ry);
            dir.rotateX(rx);

        return Ray(ori, dir.normalize());
    }
    void move(V3f m){ ori = ori + m; }
    void setRotX(float r){ rx = r; }
    void setRotY(float r){ ry = r; }
    void setRotZ(float r){ rz = r; }
    void setRot(V3f r){ rz = r.z; ry = r.y; rx = r.x; }
};

/* Light class */
enum class LightType { ambient, point, directional };
class Light{
public:
    V3f pos, dir;
    RGBA color;
    float intensity;
    LightType type;
    Light() {}
    Light(V3f pos, RGBA color, float intensity) : pos(pos), color(color), intensity(intensity), type(LightType::point) {}
    Light(RGBA color, float intensity) : color(color), intensity(intensity), type(LightType::ambient) {}
    Light(RGBA color, float intensity, V3f dir) : dir(dir), color(color), intensity(intensity), type(LightType::directional) {}
    static Light newAmbientLight(RGBA color, float intensity){ Light l; l.color = color; l.intensity = intensity; l.type = LightType::ambient; return l; }
    static Light newPointLight(V3f pos, RGBA color, float intensity){ Light l; l.pos = pos; l.color = color; l.intensity = intensity; l.type = LightType::point; return l; }
    static Light newDirLight(V3f dir, RGBA color, float intensity){ Light l; l.dir = dir; l.color = color; l.intensity = intensity; l.type = LightType::directional; return l; }
};

/* Vertex class: V3f + texture coordinates */
class Vertex {
public:
    V3f vector, normal;
    float u, v;
    Vertex() : vector(), normal(), u(0.0f), v(0.0f) {}
    Vertex(V3f vector, V3f normal, float u, float v) : vector(vector), normal(normal), u(u), v(v) {}
    void move(V3f m){ vector = vector + m; }
    void scale(float s){ vector = vector * s; normal = normal * s; }
    void scale(V3f s){ vector = V3f(s.x * vector.x, s.y * vector.y, s.z * vector.z); }
    void rotateX(float r){ vector.rotateX(r); normal.rotateX(r); }
    void rotateY(float r){ vector.rotateY(r); normal.rotateY(r); }
    void rotateZ(float r){ vector.rotateZ(r); normal.rotateZ(r); }
    void rotate(V3f r){ rotateZ(r.z); rotateY(r.y); rotateX(r.x); }
    string toString(){
        return "{" + vector.toString() + ", " + to_string(u) + ", " + to_string(v) + "}";
    }
};

/* Hit struct: hit info */
class Hit{
public:
    Vertex point;
    V3f normal, phongNormal;
    Ray ray;
    uint triId;
    Hit(){}
};

/* Material class */
class Material{
public:
    RGBA* texture, * bump;
    RGBA color;
    float diff, spec, reflective, refractive;
    Material() : texture(NULL), bump(NULL), diff(0.0f), spec(0.0f), reflective(0.0f), refractive(0.0f) {}
    Material(float diff, float spec, float reflective, float refractive) : texture(NULL), bump(NULL), diff(diff), spec(spec), reflective(reflective), refractive(refractive) {}
    void loadTexture(const char* path){ texture = loadPNG(path); }  // Load png as texture
    void loadBump(const char* path){ bump = loadPNG(path); }    // Load png as normal map
};

/* Tri class */
class Tri {
public:
    Vertex a, b, c;
    uint16_t matId;     // Material and object index for this Tri
    V3f centroid;       // Tri centroid, for BVH traversal
    uint8_t flags;      // Rendering flags

    Tri(Vertex a, Vertex b, Vertex c, uint16_t matId, uint8_t flags) : a(a), b(b), c(c), matId(matId), flags(flags), centroid((a.vector+b.vector+c.vector)*0.33333f) {}
    bool intersect(Ray ray, Hit& hit){
        V3f edge1 = b.vector - a.vector, edge2 = c.vector - a.vector, h = cross(ray.dir, edge2);

        float A = dot(edge1, h);
        if(A>-EPSILON && A<EPSILON) return false;   // Ray parallel

        V3f s = ray.ori - a.vector;
        float inv = 1.0f/A, U = inv * dot(s, h);
        if(U<0 || U>1.0f) return false;

        V3f q = cross(s, edge1);
        float V = inv * dot(ray.dir, q);
        if(V<0 || (U+V)>1.0f) return false;

        float t = inv * dot(edge2, q);
        if(t>EPSILON){  // Hit!
            hit.point.vector = ray.point(t);
            hit.point.u = (1.0f-U-V) * a.u + U * b.u + V * c.u;
            hit.point.v = (1.0f-U-V) * a.v + U * b.v + V * c.v;

            //hit.normal = interpolateNormal ? (a.normal*(1.0f-U-V) + b.normal*U + c.normal*V).normalize() : cross(a.vector - ray.point(t), b.vector - ray.point(t)).normalize();
            
            hit.phongNormal = (a.normal*(1.0f-U-V) + b.normal*U + c.normal*V).normalize();

            hit.normal = cross(a.vector - ray.point(t), b.vector - ray.point(t)).normalize();
            hit.ray = ray;
            return true;
        } else return false;
    }
    void move(V3f m){ a.move(m); b.move(m); c.move(m); }
    void scale(float s){ a.scale(s); b.scale(s); c.scale(s); }
    void scale(V3f s){ a.scale(s); b.scale(s); c.scale(s); }
    void rotate(V3f r){ a.rotate(r); b.rotate(r); c.rotate(r); }
    string toString(){
        return "{" + a.toString() + ", " + b.toString() + ", " + c.toString() + "}";
    }
};

/*  Poly class */
class Poly{
public:
    vector<Tri> tris;
    Poly() {}
    Poly(const char* path, uint matId, uint8_t flags = 0u){
        if(flags & DISABLE_RENDERING) return;
        char buff[128u];
        float x, y, z, nx, ny, nz, u, v;
        Vertex a, b, c;
        uint nvertex, nfaces;

        ifstream input; input.open(path, fstream::in);
        if(!input.is_open()){ polymsg("\e[1;91m  err opening '" + string(path) + "'\e[0m\n"); exit(EXIT_FAILURE); }

        // From the header we expect to get the number of faces and vertices, and also check if the vertex data contains the UVs and the normals
        string element;
        uint aux = 0;
        while(element!="end_header" && input.getline(buff, 128u)){
            stringstream iss(buff); iss >> element;
            if(element=="property") aux++;
            else if(element=="element"){
                iss >> element;
                if(element=="vertex"){ iss>>element; nvertex = stoi(element); }
                else { iss>>element; nfaces = stoi(element); }
            }
        }
        if(aux<8){ polymsg("\e[1;91m  err parsing '" + string(path) + "': property missing!\e[0m\n"); exit(EXIT_FAILURE); }

        // Next we get nvertex lines with the vertex data, followed by nfaces lines with the tri data
        vector<Vertex> vertices;
        for(uint i=0; i<nvertex+nfaces; i++){
            input.getline(buff, 128u);
            istringstream iss(buff);
            if(i<nvertex){  // Read vertex data: coordinates, normals and UVs, in that order
                iss >> element; x = stof(element); iss >> element; y = stof(element); iss >> element; z = stof(element);
                iss >> element; nx = stof(element); iss >> element; ny = stof(element); iss >> element; nz = stof(element);
                iss >> element; u = stof(element); iss >> element; v = stof(element);

                Vertex vertex = BLENDER ? Vertex(V3f(x,z,y), V3f(nx,nz,ny), u, 1.0f-v) : Vertex(V3f(x,y,z), V3f(nx,ny,nz), u, v);
                vertices.push_back(vertex);

            } else {        // Read tri data. Check if the line starts whith a 3
                iss >> element;
                if(element!="3") { polymsg("\e[1;91m  err parsing '" + string(path) + "': expected tri faces!\e[0m\n"); exit(EXIT_FAILURE); }
                iss >> element; a = vertices[stoi(element)];
                iss >> element; b = vertices[stoi(element)];
                iss >> element; c = vertices[stoi(element)];

                Tri tri = BLENDER ? Tri(a, c, b, matId, flags) : Tri(a, b, c, matId, flags);
                tris.push_back(tri);
            }
        }

        input.close();
        if(input.is_open()){ polymsg("\e[1;91m  err closing '" + string(path) + "'\e[0m\n"); exit(EXIT_FAILURE); }
        polymsg("\e[1;93m  loaded " + to_string(tris.size()) + " tris from '" + string(path) + "'\e[0m\n");
    }
    void move(V3f m){
        #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
        for(Tri* t=tris.data(); t<&(*tris.end()); t++) t->move(m);
    }
    void scale(V3f s){
        #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
        for(Tri* t=tris.data(); t<&(*tris.end()); t++) t->scale(s);
    }
    void scale(float s){
        #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
        for(Tri* t=tris.data(); t<&(*tris.end()); t++) t->scale(s);
    }
    void rotate(V3f r){
        #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
        for(Tri* t=tris.data(); t<&(*tris.end()); t++) t->rotate(r);
    }
    void rotateX(float r){
        #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
        for(Tri* t=tris.data(); t<&(*tris.end()); t++) t->rotate(V3f(r, 0.0f, 0.0f));
    }
    void rotateY(float r){
        #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
        for(Tri* t=tris.data(); t<&(*tris.end()); t++) t->rotate(V3f(0.0f, r, 0.0f));
    }
    void rotateZ(float r){
        #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
        for(Tri* t=tris.data(); t<&(*tris.end()); t++) t->rotate(V3f(0.0, 0.0f, r));
    }
};

/* BHV acceleration structure 
https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
*/
class BVHNode{
    V3f aabbMin, aabbMax;
    uint leftOrFirstTri, ntris; 
};

/* Renderer constructor */
PolyRenderer::PolyRenderer(){
    frame = (RGBA*) malloc(sizeof(RGBA) * WIDTH * HEIGHT);
    memset((void*) frame, 0, sizeof(RGBA*) * WIDTH * HEIGHT);
    tris = vector<Tri>(); mats = vector<Material>(); lights = vector<Light>();
}

PolyRenderer::~PolyRenderer(){
    free(frame);
}

/* Parse V3f from Node */
V3f parseV3f(YAML::Node node){
    if(node.IsSequence() && node.size()==3){
        return V3f(node[0].as<float>(), node[1].as<float>(), node[2].as<float>());
    } else throw runtime_error("\e[1;91m  err parsing V3f\e[0m\n");
}

/* Parse RGB or RGBA */
RGBA parseColor(YAML::Node node){
    if(node.IsSequence() && node.size()>2){
        return RGBA(node[0].as<uint16_t>(), node[1].as<uint16_t>(), node[2].as<uint16_t>(), node[3] ? node[3].as<uint16_t>() : 255u);
    } else throw runtime_error("\e[1;91m  err parsing RGBA\e[0m\n");
}

/* Parse rendering flags from Node */
uint8_t parseFlag(YAML::Node node){
    string flag = node.as<string>();
    if(flag=="DISABLE_RENDERING")
        return DISABLE_RENDERING;
    if(flag=="DISABLE_SHADING")
        return DISABLE_SHADING;
    if(flag=="DISABLE_TEXTURES")
        return DISABLE_TEXTURES;
    if(flag=="DISABLE_BUMP")
        return DISABLE_BUMP;
    polymsg("\e[1;96m  err unknown rendering flag '" + string(flag) + "'\e[0m\n"); return 0u;
}

/* Load scene from .poly script */
bool PolyRenderer::loadScene(const char* path){

    printf("\e[1;93m compiling '\e[95m%s\e[93m'\e[0m\n", path);

    tris.clear(); mats.clear(); lights.clear();

    string script_path = filesystem::path(path).parent_path(); script_path += "/";
    try{
        YAML::Node file = YAML::LoadFile(path);
        vector<string> mats_names, objs_names;
        vector<Poly> polys;

        // Parse camera
        if(file["camera"]){
            YAML::Node camera = file["camera"];
            V3f pos = parseV3f(camera["position"]);
            float fov = camera["fov"].as<float>();
            V3f rot = camera["rotation"] ? parseV3f(camera["rotation"]) : V3f();
            cam = new Camera(pos, fov);
            cam->setRot(rot);
        } else { polymsg("\e[1;91m  err parsing scene: camera missing\e[0m\n"); return false; }

        // Parse materials and store material name for object declaration
        if(file["materials"]){
            YAML::Node mats = file["materials"];
            this->mats.push_back(Material()); // Dummy material
            mats_names.push_back("");
            for(const auto& m : mats){
                if(m["name"] && m["diffuse"] && m["specular"]){
                    string mat_name = m["name"].as<string>();
                    float diff = m["diffuse"].as<float>(), spec = m["specular"].as<float>();
                    float reflect = m["reflective"] ? m["reflective"].as<float>() : 0.0f, refract = m["refractive"] ? m["refractive"].as<float>() : 0.0f;
                    Material mat = Material(diff, spec, reflect, refract);
                    if(m["texture"]) mat.loadTexture((script_path + m["texture"].as<string>()).c_str()); 
                    else if(m["color"]) mat.color = parseColor(m["color"]);
                    else { polymsg("\e[1;91m  err parsing material: texture or color missing!\e[0m\n"); return false; }
                    if(m["bump"]) mat.loadBump((script_path + m["bump"].as<string>()).c_str());

                    // Push material and name at same index
                    this->mats.push_back(mat);
                    mats_names.push_back(mat_name);
                } else { polymsg("\e[1;91m  err parsing material: attributes missing!\e[0m\n"); return false; }
            }
        } else { polymsg("\e[1;91m  err parsing scene: materials missing\e[0m\n"); return false; }

        // Parse objects and store object name
        uint32_t polyi = 0u;
        if(file["objects"]){
            YAML::Node objs = file["objects"];
            for(const auto& obj : objs){
                if(obj["name"] && obj["file"] && obj["material"]){
                    string obj_name = obj["name"].as<string>(), obj_file = script_path + obj["file"].as<string>(), obj_mat = obj["material"].as<string>();

                    // Get material index from name
                    auto it = find_if(mats_names.begin(), mats_names.end(),
                        [&obj_mat](const string& s){ return s==obj_mat; }
                    ); 
                    if(it==mats_names.end()){ polymsg("\e[1;91m  err parsing object: undeclared material!\e[0m\n"); return false; }
                    
                    // Parse rendering flags
                    uint8_t obj_flags = 0u;
                    YAML::Node flags = obj["flags"];
                    if(flags)
                        for(auto f : flags)
                            obj_flags |= parseFlag(f);

                    Poly poly = Poly(obj_file.c_str(), static_cast<uint16_t>(distance(mats_names.begin(), it)), obj_flags);
                    if(obj["transforms"])
                        for(const auto& t : obj["transforms"]){
                            string op = t.first.as<string>();
                            if(op=="scale"){
                                    if(t.second.IsSequence())
                                        poly.scale(parseV3f(t.second));
                                    else poly.scale(t.second.as<float>());
                                } else if(op=="rotate") poly.rotate(parseV3f(t.second));
                                else if(op=="move")
                                    poly.move(parseV3f(t.second));
                                else polymsg("\e[1;94m  unknown transform '" + op + "'\e[0m\n");
                        }
                    
                    // Insert object and object name in vector
                    polys.push_back(poly);
                    objs_names.push_back(obj_name);

                } else { polymsg("\e[1;91m  err parsing object: attributes missing!\e[0m\n"); return false; }
            }
        } else { polymsg("\e[1;91m  err parsing scene: objects missing!\e[0m\n"); return false; }

        // Once object parsing is complete, tri data is copied to the internal tri vector
        for(const auto& p : polys)
            tris.insert(tris.end(), p.tris.begin(), p.tris.end()); 

        if(file["lights"]){
            YAML::Node ls = file["lights"];
            for(const auto& l : ls){
                if(l["type"] && l["color"] && l["intensity"]){
                    string ltype = l["type"].as<string>();
                    if(ltype=="ambient"){
                        lights.push_back(Light::newAmbientLight(parseColor(l["color"]), l["intensity"].as<float>()));
                    } else if(ltype=="point"){
                        if(!l["position"]){ polymsg("\e[1;91m  err parsing point light: position missing!\e[0m\n"); return false; }
                        lights.push_back(Light::newPointLight(parseV3f(l["position"]), parseColor(l["color"]), l["intensity"].as<float>()));
                    } else if(ltype=="directional"){
                        if(!l["direction"]){ polymsg("\e[1;91m  err parsing directional light: direction missing!\e[0m\n"); return false; }
                        lights.push_back(Light::newDirLight(parseV3f(l["direction"]), parseColor(l["color"]), l["intensity"].as<float>()));
                    } else { polymsg("\e[1;91m  err parsing light: unknown light type!\e[0m\n"); return false; }
                } else { polymsg("\e[1;91m  err parsing light: attributes missing!\e[0m\n"); return false; }
            }
        } else { polymsg("\e[1;91m  err parsing scene: lights missing!\e[0m\n"); return false; }

    } catch (const YAML::ParserException& pe){ printf("\e[1;91m exception while parsing '\e[95m%s\e[93m': %s\e[0m\n", path, pe.msg.c_str()); return false; }

    // Nice printing
    system("clear");
    printIntro();
    printf("\e[1;93m compiling '\e[95m%s\e[93m' \e[92mOK\e[0m\n", path);
    for(auto s : msgvector) { printf("%s", s.c_str());}
    return true;
}

/* Texture shader: return texture element from u,v coordinates */
Fragment PolyRenderer::texture_shader(Hit& hit){
    uint matId = tris[hit.triId].matId;
    RGBA color;
    if(mats[matId].texture){
        uint16_t tx = hit.point.u * (TEXTURE_SIZE-1); tx %= (TEXTURE_SIZE-1);
        uint16_t ty = hit.point.v * (TEXTURE_SIZE-1); ty %= (TEXTURE_SIZE-1);
        color = mats[matId].texture[tx + ty*TEXTURE_SIZE];
    } else color = mats[matId].color;
    return Fragment(color);
}

/* Bump shader: update hit.normal using bump map */
V3f PolyRenderer::bump_shader(Hit& hit){
    uint matId = tris[hit.triId].matId;
    V3f normal = hit.normal;
    if(mats[matId].bump){
        uint16_t tx = hit.point.u * (TEXTURE_SIZE-1); tx %= (TEXTURE_SIZE-1);
        uint16_t ty = hit.point.v * (TEXTURE_SIZE-1); ty %= (TEXTURE_SIZE-1);
        normal =  mats[matId].bump[tx + ty*TEXTURE_SIZE].asV3f() * 2.0f - 1.0f;
    }

    return normal;
}

/* Reflection shader: compute n reflection */
Fragment PolyRenderer::reflection_shader(Ray& ray, Hit& hit, uint n){
    Material mat = mats[tris[hit.triId].matId];
    Fragment fragment = fragment_shader(hit);
    if(n<MAX_REFLECTIONS && mat.reflective>0.0f){
        V3f rdir = ray.dir - hit.normal * (dot(ray.dir, hit.normal) * 2.0f);
        Ray rray = Ray(hit.point.vector, rdir);
        Hit rhit;
        return intersection_shader(rray, rhit) ? fragment * (1.0f-mat.reflective) + reflection_shader(rray, rhit, n+1) * mat.reflective : fragment;
    } else return fragment;
}

/* Refraction shader: compute n refraction */
Fragment PolyRenderer::refraction_shader(Ray& ray, Hit& hit, uint n){
    /*Material mat = mats[tris[hit.triId].matId];
    Fragment fragment = fragment_shader(hit);
    if(n<MAX_REFRACTIONS && (mat.refractive>0.0f || fragment.a<1.0f)){
        Ray rray = Ray(hit.point.vector, ray.dir);
        Hit rhit;
        // if (fragment.a<1.0f) return intersection_shader(rray, rhit) ? refraction_shader(ray, rhit, n+1) : fragment; 
        return intersection_shader(rray, rhit) ? fragment * (1.0f-mat.refractive) + refraction_shader(ray, rhit, n+1) * mat.refractive : fragment;
    } else return fragment;*/
    return Fragment();
}

/* Fragment shader: main shading program */
Fragment PolyRenderer::fragment_shader(Hit& hit){
    Fragment out;
    Tri tri = *(&tris[hit.triId]);

    // Texture mapping step
    Fragment texture;
    if(!(tri.flags & DISABLE_TEXTURES) && tri.matId<mats.size()){
        texture = texture_shader(hit);
    } else texture = Fragment(1.0f, 1.0f, 1.0f, 1.0f);

    // TODO: bump mapping step
    if(!(tri.flags & DISABLE_BUMP) && tri.matId<mats.size()) hit.normal = bump_shader(hit);

    // Perform transparency blending
    if(texture.a<1.0f){
        Ray rray = Ray(hit.point.vector, hit.ray.dir);
        Hit rhit;
        if(intersection_shader(rray, rhit)) return texture * texture.a + fragment_shader(rhit);
    }

    // Blinn-Phong shading step
    V3f blinn_phong, ambient;
    if(!(tri.flags & DISABLE_SHADING)){
        V3f view = (cam->ori - hit.point.vector).normalize();
        for(auto& l : lights){
            if(l.intensity==0.0f) continue;

            Material mat = *(&mats[tri.matId]); 
            V3f ldir, half;
            float att;

            // Compute shading depending of light source type
            if(l.type==LightType::directional) { ldir = l.dir.normalize(); att = 1.0f; }
            else if(l.type==LightType::point) { 
                ldir = (l.pos-hit.point.vector).normalize(); 
                float dist = (l.pos-hit.point.vector).length();
                att = 1.0f / (1.0f + 0.14f * dist + 0.07 * (dist * dist));
            }
            else if(l.type==LightType::ambient) { ambient = ambient + (l.color.asV3f() * l.intensity); continue;  }

            // Check if there's geometry between the fragment and the light
            V3f sori = dot(ldir, hit.normal) < 0.0f ? hit.point.vector - hit.normal * 1e-3f : hit.point.vector + hit.normal * 1e-3f;
            Ray lray = Ray(sori, ldir);
            Hit hit2;
            if(intersection_shader(lray, hit2, (DISABLE_SHADING<<16u))) continue;
            
            // Compute diffuse and specular components of the fragment
            half = (ldir + view).normalize();
            float diff = max(0.0f, dot(hit.phongNormal, ldir));
            V3f diffuse = (l.color.asV3f() * l.intensity) * (diff * mat.diff);

            float spec = powf(max(0.0f, dot(hit.phongNormal, half)), mat.spec);
            V3f specular = (l.color.asV3f() * l.intensity) * (spec * mat.spec);

            blinn_phong = blinn_phong + diffuse*att + specular*att;
        }

    } else blinn_phong = V3f(1.0f, 1.0f, 1.0f);

    // Compute final fragment color
    out = Fragment(blinn_phong + ambient) * texture;

    return out;
}

/* Intersection shader: compute closest Tri hit */
bool PolyRenderer::intersection_shader(Ray& ray, Hit& hit, uint32_t DISCARD){
    bool intersection = false;
    float z = 1000.0f;
    Hit aux;
    for(uint i=0; i<tris.size(); i++){
        Tri tri = tris[i];
        /* DISCARD is a 4 byte flag that filters which tris should be discarted:
            - Lower 2 bytes: by material id
            - Upper first byte: unused
            - Upper second byte: by rendering flag
        */
        if(tri.matId == (uint16_t)(DISCARD & 0xffffu)) continue;
        if(tri.flags & ((DISCARD>>16u) & 0xff)) continue;

        if(tri.intersect(ray, aux) && aux.point.vector.z<z){
            hit = aux;
            z = hit.point.vector.z;
            hit.triId = i;
            intersection = true;
        }
    }
    return intersection;
}

/* Compute pixel: main pixel rendering entrypoint */
RGBA PolyRenderer::compute_pixel(uint x, uint y){
    RGBA out;
    Ray ray = cam->rayTo(x, y);
    Hit hit;
    uint reflectionRays = 0u, refractionRays = 0u;

    // Intersection step: compute closest tri intersection
    if(intersection_shader(ray, hit)){          
        Tri tri = *(&tris[hit.triId]);
        Material mat = *(&mats[tri.matId]);

        // Fragment shading step: compute fragment base color
        Fragment fragment = fragment_shader(hit);

        //return (Fragment(hit.normal, 1.0f)).toRGBA();

        // Reflection step: compute reflected color
        Fragment reflection;
        if(!(tri.flags & DISABLE_SHADING) && mat.reflective>0.0f) reflection = reflection_shader(ray, hit, 0);

        // Refraction step
        Fragment refraction;
        //if(!(tri.flags & DISABLE_SHADING) && (mat.refractive>0.0f || fragment.a<1.0f)) refraction = refraction_shader(ray, hit, 0);

        // Final pixel color
        out = (fragment + reflection + refraction).toRGBA();
    }
    
    // Return final pixel color
    return out;
}

/* Render loaded scene */
bool PolyRenderer::render(uint threads){
    float tini, trender;

    printf("\e[1;93m rendering in %s × %d\e[0m\n", getCpuCode(), threads);

    tini = static_cast<float>(omp_get_wtime());
    #pragma omp parallel for collapse(2) shared(frame, tris, mats) num_threads(threads) schedule(dynamic)
        for(int y=0; y<HEIGHT; y++)
            for(int x=0; x<WIDTH; x++)
                frame[x + y*WIDTH] = compute_pixel(x, y);
        
    trender = static_cast<float>(omp_get_wtime()) - tini;
    printf("\e[1;96m %u \e[93mtris in \e[95m%.3lfs \e[92mOK\e[m\n", static_cast<uint>(tris.size()), trender);

    return true;
}
/* Save render to .png */
void PolyRenderer::save(const char* path){
    savePNG(path, frame);
    printf("\e[1;93m render saved in '\e[95m%s\e[93m'\e[0m\n", path);
}