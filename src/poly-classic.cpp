/*
    POLY-classic ~ main header class
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include "../include/poly-classic.h"
using namespace std;

std::vector<std::string> parsermsg;

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
    printf("\e[1;93m   loaded texture '%s'\e[0m\n", path);
    parsermsg.push_back("\e[1;93m   loaded texture '" + string(path) + "'\e[0m\n");
    
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
    float fov, rx, ry, rz;
    Camera() : ori(), fov(FOV * ALPHA), rx(0.0f), ry(0.0f), rz(0.0f) {}
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
    V3f centroid;       // Centroid for BVH
    uint8_t flags;      // Rendering flags

    Tri() : a(), b(), c(), matId(0u), centroid() {}
    Tri(Vertex a, Vertex b, Vertex c) : a(a), b(b), c(c), matId(0u), centroid((a.vector+b.vector+c.vector)*0.33333f) {}
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
        char line[64], element[32];
        float x, y, z, u, v;
        Vertex a, b, c;
        int nvertex, nfaces;

        ifstream input; input.open(path, fstream::in);
        if(!input.is_open()){ printf("\e[1;91m   err opening '%s'\e[0m\n", path); exit(EXIT_FAILURE); }

        // Read vertex and face number
        for(int i=0; i<12; i++){
            input.getline(line, 64);
            stringstream iss(line);
            if(i==3){ iss >> element; iss >> element; iss >> element; nvertex = stoi(element); }
            if(i==9){ iss >> element; iss >> element; iss >> element; nfaces = stoi(element); }
        }

        // Read vertex data and store in temporal array
        Vertex* vertices = (Vertex*) malloc(sizeof(Vertex) * nvertex);
        for(uint i=0; i<nvertex; i++){
            input.getline(line, 64);
            stringstream iss(line);
            iss >> element; x = stof(element);
            iss >> element; y = stof(element);
            iss >> element; z = stof(element);
            iss >> element; u = stof(element);
            iss >> element; v = stof(element);
            
            vertices[i] = BLENDER ? Vertex(x, z, y, u, 1.0f-v) : Vertex(x, y, z, u, v);
        }

        // Create tris
        tris = vector<Tri>();
        for(uint i=0; i<nfaces; i++){
            input.getline(line, 64);
            stringstream iss(line); iss >> element;
            iss >> element; a = vertices[stoi(element)];
            iss >> element; b = vertices[stoi(element)];
            iss >> element; c = vertices[stoi(element)];

            // iss >> element; flags = individual tri rendering flags
            // iss >> element; matId = individual tri material

            Tri t = BLENDER ? Tri(a, c, b) : Tri(a, b, c); t.flags = flags; t.matId = matId;
            tris.push_back(t);
        }
        free(vertices);
        input.close();
        if(input.is_open()){ printf("\e[1;91m   err closing '%s'\e[0m\n", path); exit(EXIT_FAILURE); }

        printf("\e[1;93m   loaded %d tris from '%s'\e[0m\n", nfaces, path);
        parsermsg.push_back(("\e[1;93m   loaded " + to_string(tris.size()) + " tris from '" + string(path) + "'\e[0m\n"));
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
    tris = vector<Tri>(); mats = vector<Material>();
}

PolyRenderer::~PolyRenderer(){
    free(frame);
}

/* Parse V3f from Node */
V3f parseV3f(YAML::Node node){
    if(node.IsSequence() && node.size()==3){
        return V3f(node[0].as<float>(), node[1].as<float>(), node[2].as<float>());
    } else throw runtime_error("\e[1;91m err parsing V3f\e[0m\n");
}

/* Parse RGB or RGBA */
RGBA parseColor(YAML::Node node){
    if(node.IsSequence() && node.size()>2){
        return RGBA(node[0].as<uint16_t>(), node[1].as<uint16_t>(), node[2].as<uint16_t>(), node[3] ? node[3].as<uint16_t>() : 255u);
    } else throw runtime_error("\e[1;91m err parsing RGBA\e[0m\n");
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
    
    printf("\e[1;91m err unknown rendering flag '%s'\e[0m\n", flag.c_str()); return 0u;
}



/* Load scene from .poly script */
bool PolyRenderer::loadScene(const char* path){

    printf("\e[1;93m compiling '\e[95m%s\e[93m'\e[0m\n", path);

    tris.clear(); mats.clear();

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
            cam = new Camera(pos, fov);
        } else { printf("\e[1;91m err parsing scene: camera missing!\e[0m\n"); return false; }

        // Parse materials and store material name for object declaration
        if(file["materials"]){
            YAML::Node mats = file["materials"];
            for(const auto& m : mats){
                if(m["name"] && m["diffuse"] && m["specular"]){
                    string mat_name = m["name"].as<string>();
                    float diff = m["diffuse"].as<float>(), spec = m["specular"].as<float>();
                    Material mat = Material(diff, spec);
                    if(m["texture"]) mat.loadTexture((script_path + m["texture"].as<string>()).c_str()); 
                    else if(m["color"]) mat.color = parseColor(m["color"]);
                    else { printf("\e[1;91m err parsing material: texture or color missing!\e[0m\n"); return false; }
                    if(m["normal"]) mat.loadNormal((script_path + m["normal"].as<string>()).c_str());

                    // Push material and name at same index
                    this->mats.push_back(mat);
                    mats_names.push_back(mat_name);
                } else { printf("\e[1;91m err parsing material: attributes missing!\e[0m\n"); return false; }
            }
        } else { printf("\e[1;91m err parsing scene: materials missing!\e[0m\n"); return false; }

        // Parse objects and store object name
        if(file["objects"]){
            YAML::Node objs = file["objects"];
            for(const auto& obj : objs){
                if(obj["name"] && obj["file"] && obj["material"]){
                    string obj_name = obj["name"].as<string>(), obj_file = script_path + obj["file"].as<string>(), obj_mat = obj["material"].as<string>();

                    // Get material index from name
                    auto it = find_if(mats_names.begin(), mats_names.end(),
                        [&obj_mat](const string& s){ return s==obj_mat; }
                    ); 
                    if(it==mats_names.end()) { printf("\e[1;91m err parsing object: undeclared material!\e[0m\n"); return false; }
                    
                    // Parse rendering flags
                    uint8_t obj_flags = 0u;
                    YAML::Node flags = obj["flags"];
                    if(flags)
                        for(auto f : flags)
                            obj_flags |= parseFlag(f);

                    Poly poly = Poly(obj_file.c_str(), distance(mats_names.begin(), it), obj_flags);
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
                                else printf("\e[1;94m unknown transform '%s'\e[0m\n", op.c_str());
                        }
                    
                    // Insert object and object name in vector
                    polys.push_back(poly);
                    objs_names.push_back(obj_name);


                } else { printf("\e[1;91m err parsing object: attributes missing!\e[0m\n"); return false; }
            }
        } else { printf("\e[1;91m err parsing scene: objects missing!\e[0m\n"); return false; }

        // Once scene parsing is complete, tri data is copied to the internal tri vector
        for(const auto& p : polys)
            tris.insert(tris.end(), p.tris.begin(), p.tris.end()); 

    } catch (const YAML::ParserException& pe){ printf("\e[1;91m exception while parsing '\e[95m%s\e[93m': %s\e[0m\n", path, pe.msg.c_str()); return false; }

    // Nice printing
    system("clear");
    printIntro();
    printf("\e[1;93m compiling '\e[95m%s\e[93m' \e[92mOK\e[0m\n", path);
    for(auto s : parsermsg) { printf("%s", s.c_str());}
    return true;
}

/* Pixel shader: given a hit and a Tri, compute hit color */
RGBA PolyRenderer::pixel_shader(Hit& hit, uint triId){
    RGBA out;

    Tri tri = tris[triId];

    // Texture mapping step: returns RGBA
    if(!(tri.flags & DISABLE_TEXTURES) && tri.matId<mats.size()){
        if(mats[tri.matId].texture){
            int tx = hit.u * (TEXTURE_SIZE-1); tx %= (TEXTURE_SIZE-1);
            int ty = hit.v * (TEXTURE_SIZE-1); ty %= (TEXTURE_SIZE-1);
            out = mats[tri.matId].texture[tx + ty*TEXTURE_SIZE];
        } else out = mats[tri.matId].color;
    } else out = RGBA(255,0,255,255);   // No texture default

    // TODO: bump mapping step
    if(!(tri.flags & DISABLE_BUMP)){

    }

    // TODO: Shading step
    if(!(tri.flags & DISABLE_SHADING)){

    }
    
    return out;
}

/* Intersection shader: compute scene intersection for pixel x,y and return pixel color */
RGBA PolyRenderer::intersection_shader(uint x, uint y){
    RGBA out, ps;
    Ray ray = cam->rayTo(x, y);
    Hit hit;
    
    float z = 1000.0f;
    for(uint i=0; i<tris.size(); i++){
        Tri tri = tris[i];
        if(tri.intersect(ray, hit) && hit.vector.z<=z && !(tri.flags & DISABLE_RENDERING)){
            // Get pixel shader result
            ps = pixel_shader(hit, i);

            // Update depth only if pixel shader output is a solid color
            if(ps.a==255){ z = hit.vector.z; out = RGBA(); }
            float transparency = ps.a==255 ? 1.0f : static_cast<float>(ps.a)/255.0f;

            // Compute pixel color
            out = ps * transparency + out;
        }
    }

    return out;
}

/* Fast intersection shader: intersection shader using BVH structures */
RGBA PolyRenderer::fast_intersection_shader(uint x, uint y){
    return RGBA(0,0,0,0);
}

/* Render loaded scene */
bool PolyRenderer::render(uint threads){
    float tini, trender;

    printf("\e[1;93m device: %s × %d\e[0m\n", getCpuCode(), threads);
    printf("\e[1;93m rendering ");

    tini = static_cast<float>(omp_get_wtime());
    #pragma omp parallel for collapse(2) shared(frame, tris, mats) num_threads(threads) schedule(dynamic)
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