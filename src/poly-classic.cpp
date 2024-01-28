/*
    POLY-classic ~ main header class
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "poly-classic.h"
using namespace std;

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

/* Main vector class */
class V3f {
public:
    float x, y, z;
    V3f() : x(0.0f), y(0.0f), z(0.0f) {}
    V3f(float x, float y, float z) : x(x), y(y), z(z) {}
    V3f operator+ (V3f& a) { return V3f(x+a.x, y+a.y, z+a.z); }
    V3f operator- (V3f& a) { return V3f(x-a.x, y-a.y, z-a.z); }
    V3f operator* (V3f& a) { return V3f(x*a.x, y*a.y, z*a.z); }
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
float dot(V3f& a, V3f& b){ return a.x*b.x + a.y*b.y + a.z*b.z; }  
V3f cross(V3f& a, V3f& b){ return V3f(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }  

/* RGBA pixel class */
class RGBA{
public:
    uint16_t r, g, b, a;
    RGBA() : r(0u), g(0u), b(0u), a(0u) {}
    RGBA(uint16_t r, uint16_t g, uint16_t b, uint16_t a) : r(r), g(g), b(b), a(a) {}
    //RGBA(float r, float g, float b, float a) : r(static_cast<uint16_t>(r)), g(static_cast<uint16_t>(g)), b(static_cast<uint16_t>(b)), a(static_cast<uint16_t>(a)) {}
    void clamp255(){ r = r>255 ? 255 : r; g = g>255 ? 255 : g; b = b>255 ? 255 : b; a = a>255 ? 255 : a; }
    RGBA operator+ (RGBA& color) { RGBA out = RGBA(r+color.r, g+color.g, b+color.b, a+color.a); out.clamp255(); return out; }
    RGBA operator- (RGBA& color) { RGBA out = RGBA(r-color.r, g-color.g, b-color.b, a-color.a); out.clamp255(); return out; }
    RGBA operator* (float f) { RGBA out = RGBA(static_cast<uint16_t>(r*f), static_cast<uint16_t>(g*f), static_cast<uint16_t>(b*f), static_cast<uint16_t>(a*f)); out.clamp255(); return out; }
    RGBA operator* (V3f v) { RGBA out = RGBA(static_cast<uint16_t>(r*v.x), static_cast<uint16_t>(g*v.y), static_cast<uint16_t>(b*v.z), a); out.clamp255(); return out; }
};

/* Load png as RGBA texture */
RGBA* loadPNG(const char* path){
    FILE* input = fopen(path, "rb");
        if(!input){ printf("\e[1;91m  err loading '%s': file could not be opened!\n\e[0m", path); exit(EXIT_FAILURE); }

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if(!png){ 
            fclose(input);
            std::cerr << "\e[1;91m  ERR! texture '" << path << "' read struct could not be loaded\e[0m\n"; exit(EXIT_FAILURE); 
        }

    png_infop info = png_create_info_struct(png);
        if(!info){
            fclose(input);
            png_destroy_read_struct(&png, NULL, NULL);
            std::cerr << "\e[1;91m  ERR! texture '" << path << "' info struct could not be loaded\e[0m\n"; exit(EXIT_FAILURE);
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

    if(fclose(input)==-1){ printf("\e[1;91m  err closing '%s'!\n\e[0m", path); exit(EXIT_FAILURE); }
    printf("\e[1;93m  loaded texture '%s'\e[0m\n", path);
    
    return texture;
}

/* Save RGBA* to png */
void savePNG(const char* path, RGBA* texture){
    FILE* output = fopen(path, "wb");
        if(!output){ printf("\e[1;91m err saving to '%s': file could not be opened!\e[0m\n", path); exit(EXIT_FAILURE); }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if(!png){ 
            fclose(output);
            cerr << "\e[1;91mERR! error saving render\e[0m\n"; exit(EXIT_FAILURE);
        }

    png_infop info = png_create_info_struct(png);
        if(!info){ 
            fclose(output);
            png_destroy_write_struct(&png, NULL);
            cerr << "\e[1;91mERR! error saving render\e[0m\n"; exit(EXIT_FAILURE);
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

    if(fclose(output)==-1){ printf("\e[1;91m err closing '%s'\e[0m\n", path); exit(EXIT_FAILURE); }
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
    float ar, fov, rx, ry, rz;
    Camera() : ori(), ar(AR), fov(FOV * ALPHA), rx(0.0f), ry(0.0f), rz(0.0f) {}
    Camera(V3f ori, float ar, float fov) : ori(ori), ar(ar), fov(fov*ALPHA), rx(0.0f), ry(0.0f), rz(0.0f) {}
    Ray rayTo(int x, int y){
        float aux = tanf(fov / 2.0f);
        float px = (2.0f * ((x + 0.5f)/WIDTH) - 1.0f) * aux * ar;
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
    void setRot(V3f r){ rz = r.z; ry = ry; rx = r.x; }
};

/* Vertex class: V3f + texture coordinates */
class Vertex {
public:
    V3f vector;
    float u, v;
    Vertex() : vector(), u(0.0f), v(0.0f) {}
    Vertex(float x, float y, float z, float u, float v) : vector(x, y, z), u(u), v(v) {}
    void move(V3f m){ vector = vector + m; }
    void scale(float s){ vector = vector * s; }
    void scale(V3f s){ vector = V3f(s.x * vector.x, s.y * vector.y, s.z * vector.z); }
    void rotateX(float r){ vector.rotateX(r); }
    void rotateY(float r){ vector.rotateY(r); }
    void rotateZ(float r){ vector.rotateZ(r); }
    void rotate(V3f r){ rotateZ(r.z); rotateY(r.y); rotateX(r.x); }
    string toString(){
        return "{" + vector.toString() + ", " + to_string(u) + ", " + to_string(v) + "}";
    }
};

typedef Vertex Hit;

/* Material class */
class Material{
public:
    RGBA* texture, * normal;
    RGBA color;
    float diff, spec;
    Material() : texture(NULL), normal(NULL), diff(0.0f), spec(0.0f) {}
    Material(float diff, float spec) : texture(NULL), normal(NULL), diff(diff), spec(spec) {}
    void loadTexture(const char* path){ texture = loadPNG(path); }  // Load png as texture
    void loadNormal(const char* path){ normal = loadPNG(path); }    // Load png as normal map
};

/* Tri class */
class Tri {
public:
    Vertex a, b, c;
    uint matId;         // Material index for this Tri
    Tri() : a(), b(), c() {}
    Tri(Vertex a, Vertex b, Vertex c) : a(a), b(b), c(c) {}
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
            hit.vector = ray.point(t);
            hit.u = (1.0f-U-V) * a.u + U * b.u + V * c.u;
            hit.v = (1.0f-U-V) * a.v + U * b.v + V * c.v;
            return true;
        } else return false;
    }
    void move(V3f m){ a.move(m); b.move(m); c.move(m); }
    void scale(float s){ a.scale(s); b.scale(s); c.scale(s); }
    void scale(V3f s){ a.scale(s); b.scale(s); c.scale(s); }
    void rotate(V3f r){ a.rotate(r); b.rotate(r); c.rotate(r); }
};

/* Renderer constructor */
PolyRenderer::PolyRenderer(){
    frame = (RGBA*) malloc(sizeof(RGBA) * WIDTH * HEIGHT);
    memset((void*) frame, 0, sizeof(RGBA*) * WIDTH * HEIGHT);
    tris = 0; mats = 0;
}

PolyRenderer::~PolyRenderer(){
    free(frame);
    if(tris) free(tris);
    if(mats) free(mats);
}

/* Load scene from .poly script */
bool PolyRenderer::loadScene(const char* path){

    vector<Tri> trivec;
    vector<Material> matvec;

    cam = new Camera(V3f(0.0f, 0.0f, 0.0f), 4.0f/3.0f, 90.0f);

    Tri t = Tri(Vertex(0.0f, 0.0f, 0.0f, 0.0f, 1.0f), Vertex(0.0f, 1.0f, 0.0f, 0.0f, 0.0f), Vertex(1.0f, 0.0f, 0.0f, 1.0f, 1.0f));
    t.matId = 0;
    t.move(V3f(-0.5f, 0.0f, 1.2f));
    trivec.push_back(t);

    t = Tri(Vertex(0.0f, 1.0f, 0.0f, 0.0f, 0.0f), Vertex(1.0f, 1.0f, 0.0f, 1.0f, 0.0f), Vertex(1.0f, 0.0f, 0.0f, 1.0f, 1.0f));
    t.matId = 0;
    t.move(V3f(-0.5f, 0.0f, 1.2f));
    trivec.push_back(t);

    Material m = Material(2.0f, 32.0f);
    m.loadTexture("demo/beso.png");
    matvec.push_back(m);

    // Copy data to renderer
    N_TRIS = trivec.size(); N_MATS = matvec.size();
    tris = (Tri*) malloc(sizeof(Tri)*N_TRIS); copy(trivec.begin(), trivec.end(), tris);
    mats = (Material*) malloc(sizeof(Material)*N_MATS); copy(matvec.begin(), matvec.end(), mats);

    printf("\e[1;93m loading '\e[95m%s\e[93m' \e[92mOK\e[0m\n", path);

    return true;
}

/* Pixel shader: given a hit and a Tri, compute pixel color */
RGBA PolyRenderer::pixel_shader(Hit& hit, uint triId){
    RGBA out;

    int tx = hit.u * (TEXTURE_SIZE-1); tx %= (TEXTURE_SIZE-1);
    int ty = hit.v * (TEXTURE_SIZE-1); ty %= (TEXTURE_SIZE-1);

    Tri tri = tris[triId];
    
    out = mats[tri.matId].texture[tx + ty*TEXTURE_SIZE];

    return out;
}

/* Intersection shader: compute scene intersection for pixel x,y */
RGBA PolyRenderer::intersection_shader(int x, int y){
    RGBA out;
    Ray ray = cam->rayTo(x, y);
    Hit hit;
    
    float z = 1000.0f;
    for(uint i=0; i<N_TRIS; i++){
        Tri tri = tris[i];
        if(tri.intersect(ray, hit) && hit.vector.z<=z){
            z = hit.vector.z;
            out = pixel_shader(hit, i);
        }
    }

    return out;
}

/* Render loaded scene */
bool PolyRenderer::render(uint threads){
    float tini, trender;

    printf("\e[1;93m device: %s × %d\e[0m\n", getCpuCode(), threads);
    printf("\e[1;93m rendering ");

    tini = static_cast<float>(omp_get_wtime());
    #pragma omp parallel for collapse(2) shared(frame) num_threads(threads) schedule(static) shared(tris, mats)
        for(int y=0; y<HEIGHT; y++)
            for(int x=0; x<WIDTH; x++)
                frame[x + y*WIDTH] = intersection_shader(x, y);
        
    trender = static_cast<float>(omp_get_wtime()) - tini;
    printf("\e[95m%.3lfs \e[92mOK\e[m\n", trender);

    return true;
}
/* Save render to .png */
void PolyRenderer::save(const char* path){
    savePNG(path, frame);
    printf("\e[1;93m render saved in '\e[95m%s\e[93m'\e[0m\n", path);
}