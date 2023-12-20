/*
    POLY classc ~ main header class
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
    void clamp255(){ r = r>255 ? 255 : r; g = g>255 ? 255 : g; b = b>255 ? 255 : b; a = a>255 ? 255 : a; }
    RGBA operator+ (RGBA& color) { RGBA out = RGBA(r+color.r, g+color.g, b+color.b, a+color.a); out.clamp255(); return out; }
    RGBA operator- (RGBA& color) { RGBA out = RGBA(r-color.r, g-color.g, b-color.b, a-color.a); out.clamp255(); return out; }
    RGBA operator* (float f) { RGBA out = RGBA(static_cast<uint16_t>(r*f), static_cast<uint16_t>(g*f), static_cast<uint16_t>(b*f), static_cast<uint16_t>(a*f)); out.clamp255(); return out; }
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
    void rotateX(float r){ rx = r; }
    void rotateY(float r){ ry = r; }
    void rotateZ(float r){ rz = r; }
};

/* Renderer constructor */
PolyRenderer::PolyRenderer(){
    frame = (RGBA*) malloc(sizeof(RGBA) * WIDTH * HEIGHT);
    memset((void*) frame, 0, sizeof(RGBA*) * WIDTH * HEIGHT);
}

PolyRenderer::~PolyRenderer(){
    free(frame);
}

/* Load scene from .poly script */
bool PolyRenderer::loadScene(const char* path){
    return true;
}

/* Main shader function; each thread run this when rendering pixel (x,y) */
RGBA shader(int x, int y){
    RGBA out;


    return out;
}

/* Render loaded scene */
bool PolyRenderer::render(uint threads){
    float tini, trender;

    printf("\e[1;93m device: %s × %d\e[0m\n", getCpuCode(), threads);
    printf("\e[1;93m rendering ");

    tini = static_cast<float>(omp_get_wtime());
    #pragma omp parallel for collapse(2) shared(frame) num_threads(threads) schedule(static)
        for(int y=0; y<HEIGHT; y++)
            for(int x=0; x<WIDTH; x++)
                frame[x + y*WIDTH] = shader(x, y);
        
    trender = static_cast<float>(omp_get_wtime()) - tini;
    printf("\e[95m%.3lfs \e[92mOK\e[m\n", trender);

    return true;
}
/* Save render to .png */
void PolyRenderer::save(const char* path){
    savePNG(path, frame);
    printf("\e[1;93m render saved in '\e[95m%s\e[93m'\e[0m\n", path);
}