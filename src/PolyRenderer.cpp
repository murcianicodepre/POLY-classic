#include "PolyRenderer.h"
#include "Material.h"

/*
    PolyRenderer ~ main poly-classic renderer class 
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

std::vector<std::string> msgvector;

// BVHNode class for acceleration structure
class BVHNode{
public:
    Vec3 aabbMin, aabbMax;
    uint32_t left, first, n;
    BVHNode();
    bool intersectAABB(Ray&);
};
BVHNode::BVHNode() : aabbMin(), aabbMax(), left(0), first(0), n(0) {}
bool BVHNode::intersectAABB(Ray& ray){
    // Slab test for volume-ray intersection
    float tx1 = (aabbMin.x - ray.ori.x) / ray.dir.x, tx2 = (aabbMax.x - ray.ori.x) / ray.dir.x;
    float tmin = fmin(tx1, tx2), tmax = fmax(tx1, tx2);
    float ty1 = (aabbMin.y - ray.ori.y) / ray.dir.y, ty2 = (aabbMax.y - ray.ori.y) / ray.dir.y;
    tmin = fmax(fmin(ty1, ty2), tmin), tmax = fmin(fmax(ty1, ty2), tmax);
    float tz1 = (aabbMin.z - ray.ori.z) / ray.dir.z, tz2 = (aabbMax.z - ray.ori.z) / ray.dir.z;
    tmin = fmax(fmin(tz1, tz2), tmin), tmax = fmin(fmax(tz1, tz2), tmax);
    return tmax >= tmin && tmax > 0.0f;
}

void PolyRenderer::buildBVH(){
    bvh = (BVHNode*) malloc(sizeof(BVHNode) * 2 * tris.size() - 1);
    triIdx = (uint32_t*) malloc(sizeof(uint32_t) * tris.size());
    for(uint32_t i=0; i<tris.size(); i++){
        Tri& tri = tris[i];
        triIdx[i] = i;
        tri.centroid = (tri.a.xyz + tri.b.xyz + tri.c.xyz) * 0.33333f;
    }

    // Assign all tris to root node
    BVHNode& root = bvh[0];
    root.n = tris.size(); root.first = 0; root.left = 0;

    updateNodeBounds(0);    // Update bounds of root node
    subdivide(0);           // Start subdivision

    // Print size of bvh structure
    PolyRenderer::polyMsg("\e[1;93m compiling bvh: \e[95m" + to_string(nextNode) + " \e[93mnodes \e[92mOK\e[0m\n");
}

void PolyRenderer::updateNodeBounds(uint32_t nodeId){
    BVHNode& node = bvh[nodeId];
    node.aabbMin = Vec3(1e30f), node.aabbMax = Vec3(1e-30f);
    for(uint32_t first = node.first, i=0; i<node.n; i++){
        Tri& tri = tris[triIdx[first+i]];
        node.aabbMin = Vec3::minVec3(node.aabbMin, tri.a.xyz);
        node.aabbMin = Vec3::minVec3(node.aabbMin, tri.b.xyz);
        node.aabbMin = Vec3::minVec3(node.aabbMin, tri.c.xyz);
        node.aabbMax = Vec3::maxVec3(node.aabbMax, tri.a.xyz);
        node.aabbMax = Vec3::maxVec3(node.aabbMax, tri.b.xyz);
        node.aabbMax = Vec3::maxVec3(node.aabbMax, tri.c.xyz);
    }
}

void PolyRenderer::subdivide(uint32_t nodeId){
    BVHNode& node = bvh[nodeId];
    if(node.n<=2) return;  // Terminate recursive subdivision

    // Get split axis
    Vec3 ext = node.aabbMax - node.aabbMin;
    uint8_t axis = 0u;
    if(ext.y > ext.x) axis = 1;
    if(ext.z > ext[axis]) axis = 2;
    float split = node.aabbMin[axis] + ext[axis] * 0.5f;

    // Split the geometry in two parts
    uint32_t i = node.first, j = i + node.n - 1;
    while(i<=j){
        Tri& tri = tris[triIdx[i]];
        if(tri.centroid[axis]>=split) i++;
        else std::swap(triIdx[i], triIdx[j--]);
    }

    // Terminate if one of the sides is empty
    uint32_t leftCount = i - node.first;
    if(leftCount==0 || leftCount==node.n) return;

    // Create child nodes
    uint32_t leftIdx = nextNode++, rightIdx = nextNode++;
    bvh[leftIdx].first = node.first; bvh[leftIdx].n = leftCount;
    bvh[rightIdx].first = i; bvh[rightIdx].n = node.n - leftCount;
    node.left = leftIdx; node.n = 0;
    updateNodeBounds(leftIdx); updateNodeBounds(rightIdx);

    // Continue recursion
    subdivide(leftIdx); subdivide(rightIdx);
}

void PolyRenderer::intersectBVH(Ray& ray, Hit& hit, uint32_t nodeId, bool& bvhIntersection, uint16_t flags){
    BVHNode& node = bvh[nodeId];
    if(!node.intersectAABB(ray)) return;
    bvhIntersection = true;
    if(node.n>0){
        float z = 1000.0f;
        Hit aux;
        for(uint32_t i=0; i<node.n; i++){
            Tri& tri = tris[triIdx[node.first+i]];

            // Discard tri if matches with low byte of flags
            if(tri.flags & static_cast<uint8_t>(flags & 0xffu)) continue;

            // Tri intersection test
            if(tri.intersect(ray, aux) && aux.point.z<z){
                hit.ray = ray;
                hit = aux;
                z = hit.point.z;
                hit.tri = triIdx[node.first+i];
                hit.valid = true;
            }
        }

    } else {
        intersectBVH(ray, hit, node.left, bvhIntersection, flags);
        intersectBVH(ray, hit, node.left + 1, bvhIntersection, flags);
    }
}

// PolyRenderer: main poly-classic renderer class
PolyRenderer::PolyRenderer(){
    frame = (RGBA*) malloc(sizeof(RGBA)*WIDTH*HEIGHT);
    memset((void*) frame, 0, sizeof(RGBA) * WIDTH * HEIGHT);
    tris = vector<Tri>(); mats = vector<Material>(); lights = vector<unique_ptr<Light>>();
}
PolyRenderer::~PolyRenderer(){
    free(frame);
    free(bvh);
    free(triIdx);
}

// Scene parser functions
Vec3 PolyRenderer::parseVec3(YAML::Node node){
    if(node.IsSequence() && node.size()==3){
        return Vec3(node[0].as<float>(), node[1].as<float>(), node[2].as<float>());
    } else throw runtime_error("\e[1;91m  err parsing Vec3\e[0m\n");
}
RGBA PolyRenderer::parseColor(YAML::Node node){
    if(node.IsSequence() && node.size()>2){
        return RGBA(node[0].as<uint8_t>(), node[1].as<uint8_t>(), node[2].as<uint8_t>(), node[3] ? node[3].as<uint8_t>() : 255u);
    } else throw runtime_error("\e[1;91m  err parsing RGBA\e[0m\n");
}

uint16_t PolyRenderer::parseFlags(YAML::Node node){
    uint16_t flags = 0x0000u;
    for(auto f : node){
        string flag = f.as<string>();
        if(flag=="DISABLE_RENDERING")
            flags |= DISABLE_RENDERING;
        else if(flag=="DISABLE_SHADING")
            flags |= DISABLE_SHADING;
        else if(flag=="DISABLE_TEXTURES")
            flags |= DISABLE_TEXTURES;
        else if(flag=="DISABLE_BUMP")
            flags |= DISABLE_BUMP;
        else if(flag=="DISABLE_TRANSPARENCY")
            flags |= DISABLE_TRANSPARENCY;
        else if(flag=="DISABLE_SHADOWS")
            flags |= DISABLE_SHADOWS;
        else if(flag=="DISABLE_REFLECTIONS")
            flags |= DISABLE_REFLECTIONS;
        else if(flag=="DISABLE_REFRACTIONS")
            flags |= DISABLE_REFRACTIONS;
        else if(flag=="DRAW_AABB")
            flags |= (DRAW_AABB<<8);
        else if(flag=="DISABLE_FAST_INTERSECTION_SHADER")
            flags |= (DISABLE_FAST_INTERSECTION_SHADER<<8);
        else PolyRenderer::polyMsg("\e[1;96m  err unknown flag '" + string(flag) + "'\e[0m\n");
        
    }
    return flags;
}

// Load scene data from .poly
bool PolyRenderer::loadScene(const char* path){
    
    printf("\e[1;93m compiling '\e[95m%s\e[93m'\e[0m\n", path);

    tris.clear(); mats.clear(); lights.clear();

    string script_path = filesystem::path(path).parent_path(); if(!script_path.empty()) script_path += "/";
    try{
        YAML::Node file = YAML::LoadFile(path);
        vector<string> mats_names, objs_names;
        vector<Poly> polys;

        // Parse debug flags
        if(file["debug"])
            debug = (parseFlags(file["debug"]) & 0xffffu);
        
        // Parse camera
        if(file["camera"]){
            YAML::Node camera = file["camera"];
            Vec3 pos = parseVec3(camera["position"]);
            float fov = camera["fov"].as<float>();
            Vec3 rot = camera["rotation"] ? parseVec3(camera["rotation"]) : Vec3();
            cam = new Camera(pos, rot, fov);
        } else { PolyRenderer::polyMsg("\e[1;91m  err parsing scene: camera missing\e[0m\n"); return false; }

        // Parse materials and store material name for object declaration
        if(file["materials"]){
            YAML::Node mats = file["materials"];
            for(const auto& m : mats){
                if(m["name"] && m["diffuse"] && m["specular"]){
                    string mat_name = m["name"].as<string>();
                    float diff = m["diffuse"].as<float>(), spec = m["specular"].as<float>();
                    float reflect = m["reflective"] ? m["reflective"].as<float>() : 0.0f, refract = m["refractive"] ? m["refractive"].as<float>() : 0.0f;
                    Material mat = Material(diff, spec, reflect, refract);
                    if(m["texture"]) mat.loadTexture((script_path + m["texture"].as<string>()).c_str()); 
                    else if(m["color"]) mat.color = parseColor(m["color"]);
                    else { PolyRenderer::polyMsg("\e[1;91m  err parsing material: texture or color missing!\e[0m\n"); return false; }
                    if(m["bump"]) mat.loadBump((script_path + m["bump"].as<string>()).c_str());

                    // Push material and name at same index
                    this->mats.push_back(mat);
                    mats_names.push_back(mat_name);
                } else { PolyRenderer::polyMsg("\e[1;91m  err parsing material: attributes missing!\e[0m\n"); return false; }
            }
        } else { PolyRenderer::polyMsg("\e[1;91m  err parsing scene: materials missing\e[0m\n"); return false; }

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
                    if(it==mats_names.end()){ PolyRenderer::polyMsg("\e[1;91m  err parsing object: undeclared material!\e[0m\n"); return false; }
                    
                    // Parse rendering flags
                    uint8_t obj_flags = 0u;
                    if(obj["flags"])
                        obj_flags = static_cast<uint8_t>(parseFlags(obj["flags"]));

                    Poly poly = Poly(obj_file.c_str(), static_cast<uint16_t>(distance(mats_names.begin(), it)), (obj_flags|(debug & 0xffu)));
                    if(obj["transforms"])
                        for(const auto& t : obj["transforms"]){
                            string op = t.first.as<string>();
                            if(op=="scale"){
                                    if(t.second.IsSequence())
                                        poly.scale(parseVec3(t.second));
                                    else poly.scale(t.second.as<float>());
                                } else if(op=="rotate") poly.rotate(parseVec3(t.second));
                                else if(op=="move")
                                    poly.move(parseVec3(t.second));
                                else PolyRenderer::polyMsg("\e[1;94m  unknown transform '" + op + "'\e[0m\n");
                        }

                    // Insert both object and object name in vector
                    polys.push_back(poly);
                    objs_names.push_back(obj_name);

                } else { PolyRenderer::polyMsg("\e[1;91m  err parsing object: attributes missing!\e[0m\n"); return false; }
            }
        } else { PolyRenderer::polyMsg("\e[1;91m  err parsing scene: objects missing!\e[0m\n"); return false; }

        // Once object parsing is complete, tri data is copied to the internal tri vector
        for(const auto& p : polys)
            tris.insert(tris.end(), p.tris.begin(), p.tris.end()); 

        if(file["lights"]){
            YAML::Node ls = file["lights"];
            for(const auto& l : ls){
                if(l["type"] && l["color"] && l["intensity"]){
                    string ltype = l["type"].as<string>();
                    if(ltype=="ambient"){
                        lights.push_back(make_unique<Ambient>(parseColor(l["color"]), l["intensity"].as<float>()));
                    } else if(ltype=="point"){
                        if(!l["position"]){ PolyRenderer::polyMsg("\e[1;91m  err parsing point light: position missing!\e[0m\n"); return false; }
                        lights.push_back(make_unique<PointLight>(parseVec3(l["position"]), parseColor(l["color"]), l["intensity"].as<float>()));
                    } else if(ltype=="directional"){
                        if(!l["direction"]){ PolyRenderer::polyMsg("\e[1;91m  err parsing directional light: direction missing!\e[0m\n"); return false; }
                        lights.push_back(make_unique<DirectionalLight>(parseVec3(l["direction"]), parseColor(l["color"]), l["intensity"].as<float>()));
                    } else { PolyRenderer::polyMsg("\e[1;91m  err parsing light: unknown light type!\e[0m\n"); return false; }
                } else { PolyRenderer::polyMsg("\e[1;91m  err parsing light: attributes missing!\e[0m\n"); return false; }
            }
        } else { PolyRenderer::polyMsg("\e[1;91m  err parsing scene: lights missing!\e[0m\n"); return false; }

    } catch (const YAML::ParserException& pe){ printf("\e[1;91m exception while parsing '\e[95m%s\e[93m': %s\e[0m\n", path, pe.msg.c_str()); return false; }

    // Nice printing
    system("clear");
    printIntro();
    printf("\e[1;93m compiling '\e[95m%s\e[93m' \e[92mOK\e[0m\n", path);
    for(auto s : msgvector) { printf("%s", s.c_str());}
    return true;
}

// Render scene using N threads
bool PolyRenderer::render(uint8_t threads){

    // If global option DISABLE_RENDERING is set, return false directly
    if((debug & 0xffu) & DISABLE_RENDERING){
        printf("\e[1;91m fatal: rendering disabled!\e[0m\n");
        return false;
    }

    float tini, trender;

    // Open display and create X11 window to show rendering at real time
    if(!XInitThreads()){ printf("\e[1;91m X11 err!\e[0m\n"); exit(EXIT_FAILURE); }
    Display* display = XOpenDisplay(nullptr);
    if(!display){ printf("\e[1;91m err opening display!\e[0m\n"); exit(EXIT_FAILURE); }
    int screen = DefaultScreen(display);
    Window win = XCreateSimpleWindow(display, RootWindow(display, screen), 0,0, WIDTH, HEIGHT, 0, BlackPixel(display, screen), WhitePixel(display, screen));
    XSelectInput(display, win, KeyPressMask | KeyReleaseMask | ButtonPressMask | ButtonReleaseMask | ExposureMask);
    XStoreName(display, win, "POLY-classic ~ rendering");
    
    // Set fixed window size
    XSizeHints sizeHints; sizeHints.flags = PMinSize | PMaxSize;
    sizeHints.min_width = WIDTH;
    sizeHints.max_width = WIDTH;
    sizeHints.min_height = HEIGHT;
    sizeHints.max_height = HEIGHT;
    XSetWMNormalHints(display, win, &sizeHints);

    XSetWindowBackground(display, win, 0xff << 16 | 0x00 << 8 | 0xff);
    XClearWindow(display, win);

    XMapWindow(display, win);
    GC gc = XCreateGC(display, win, 0, nullptr);

    Atom wmdel = XInternAtom(display, "WM_DELETE_WINDOW", true);
    XSetWMProtocols(display, win, &wmdel, 1);

    // Wait until the rendering window pops out
    XEvent event;
    while(true){
        XNextEvent(display, &event);
        if(event.type==Expose) break;
    }

    // Build BVH acceleration struct after scene is loaded
    if(!((debug>>8) & DISABLE_FAST_INTERSECTION_SHADER))
        buildBVH();

    // Start timer and launch rendering threads
    printf("\e[1;93m rendering in %s × %d\e[0m\n", PolyRenderer::getCpu(), threads);
    tini = static_cast<float>(omp_get_wtime());

    #pragma omp parallel for collapse(2) shared(frame, tris, mats, lights) num_threads(threads) schedule(dynamic)
    for(uint by=0; by<HEIGHT/TILE_SIZE; by++)
        for(uint bx=0; bx<WIDTH/TILE_SIZE; bx++){
            RGBA* tile = (RGBA*) malloc(sizeof(RGBA) * TILE_SIZE * TILE_SIZE);

            // Compute tile
            for(uint ty=0; ty<TILE_SIZE; ty++)
                for(uint tx=0; tx<TILE_SIZE; tx++){
                    tile[tx + ty*TILE_SIZE] = compute_pixel(tx + bx*TILE_SIZE, ty + by*TILE_SIZE);
                    #pragma omp critical
                    frame[(tx+bx*TILE_SIZE) + (ty+by*TILE_SIZE)*WIDTH] = tile[tx + ty*TILE_SIZE];
                }

            // Draw tile in rendering window and poll window events
            #pragma omp critical
            {
                XImage* tileImg = XCreateImage(display, DefaultVisual(display, screen), DefaultDepth(display, screen), 
                                              ZPixmap, 0, nullptr, TILE_SIZE, TILE_SIZE, 32,0);
                tileImg->data = reinterpret_cast<char*>(tile);

                XPutImage(display, win, gc, tileImg, 0,0, bx*TILE_SIZE, by*TILE_SIZE, TILE_SIZE, TILE_SIZE);
            }
            free(tile);
        }
        
    // Compute elapsed time
    trender = static_cast<float>(omp_get_wtime()) - tini;
    printf("\e[1;96m %u \e[93mtris in \e[95m%.3lfs \e[92mOK\e[m\n", static_cast<uint>(tris.size()), trender);

    XStoreName(display, win, string("POLY-classic ~ finished in " + to_string(trender) + "s").c_str());

    // Loop until window is closed
    while(true){
        XNextEvent(display, &event);
        if(event.type==KeyPress){
            KeySym key = XLookupKeysym(&event.xkey, 0);
            char keyChar = XKeysymToKeycode(display, key);
            if(keyChar==0x1b || keyChar==0x03){ XFreeGC(display, gc); XCloseDisplay(display); return true; }
        }
        else if(event.type==ClientMessage && event.xclient.data.l[0]==wmdel) { XFreeGC(display, gc); XCloseDisplay(display); return true; }
    }

    return true;
}

// Compute pixel for that (x,y)
RGBA PolyRenderer::compute_pixel(uint16_t x, uint16_t y){
    RGBA out;
    Ray ray = cam->rayTo(x, y);
    Hit hit = {};
    Fragment frag;

    // Intersection step: compute closest tri intersection
    if(intersection_shader(ray, hit)){
        
        // Debug: draw AABB
        if(((debug>>8) & DRAW_AABB) && !hit.valid) return RGBA(0,255,255,255);

        Tri& tri = tris[hit.tri];
        Material& mat = mats[tri.mat];

        // Fragment shader step: compute fragment color
        frag = fragment_shader(hit);

        // Reflection step: compute N reflection
        frag = (!(tri.flags & (DISABLE_SHADING | DISABLE_REFLECTIONS)) && mat.reflective>1e-3f) ? reflection_shader(ray, hit, 0) : frag;

        // TODO Refraction step: compute N refraction
        frag = (!(tri.flags & (DISABLE_SHADING | DISABLE_REFRACTIONS )) && mat.refractive>1e-3f) ? refraction_shader(ray, hit, 0) : frag;

        // Compute final pixel color
        out = RGBA(frag);
    }

    // Return final pixel color;
    return out;
}

// Intersection shader: compute closest tri hit
bool PolyRenderer::intersection_shader(Ray& ray, Hit& hit, uint16_t flags){

    if ((debug>>8) & DISABLE_FAST_INTERSECTION_SHADER){
        float z = 1000.0f;
        Hit aux;
        for(uint32_t i=0; i<tris.size(); i++){
            Tri& tri = tris[i];

            // Discard if tri.flags matches with low byte of flags
            if(tri.flags & static_cast<uint8_t>(flags & 0xffu)) continue;

            // Intersection test
            if(tri.intersect(ray, aux) && aux.point.z<z){
                hit = aux; z = hit.point.z; hit.tri = i; hit.valid = true; hit.ray = ray;
            }
        }
        return hit.valid;
    } else {
        bool bvhIntersection = false;
        intersectBVH(ray, hit, 0, bvhIntersection, flags);
        return ((flags>>8) & DRAW_AABB) ? bvhIntersection : hit.valid;
    }
}

// Fragment shader: main fragment shader program
Fragment PolyRenderer::fragment_shader(Hit& hit){
    Fragment out;
    Tri& tri = tris[hit.tri];
    Material& mat = mats[tri.mat];

    // Texture mapping step
    Fragment tex = (!((tri.flags | (debug & 0xffu)) & DISABLE_TEXTURES) && tri.mat<mats.size()) ? texture_shader(hit) : Fragment(1.0f,0.0f,1.0f,1.0f);

    // Bump mapping step
    Vec3 bump = (!((tri.flags | (debug & 0xffu)) & DISABLE_BUMP) && tri.mat<mats.size()) ? bump_shader(hit) : Vec3(1.0f,1.0f,1.0f);
    hit.phong = hit.phong * bump;

    // Transparency blending
    if(!((tri.flags | (debug & 0xffu)) & DISABLE_TRANSPARENCY) && tex.a<1.0f){
        Ray rray = Ray(hit.point, hit.ray.dir);
        Hit rhit;
        if(intersection_shader(rray, rhit)) return (tex * tex.a) + fragment_shader(rhit);
    }

    // Blinn-Phong shading step
    Vec3 blinn_phong;
    if(!((tri.flags | (debug & 0xffu)) & DISABLE_SHADING)){
        Vec3 view = (cam->ori - hit.point).normalize();

        for(auto& l : lights){
            if(l.get()->intensity<1e-3f) continue;
            Vec3 ldir, half;
            float att;

            // Compute ldir and att depending on the light type
            if(l.get()->type()==LightType::Point){
                PointLight* pl = dynamic_cast<PointLight*>(l.get());
                ldir = (pl->pos - hit.point).normalize();
                float dist = (pl->pos - hit.point).length();
                att = 1.0f / (1.0f + 0.14f * dist + 0.07f * (dist*dist));
            } else if(l.get()->type()==LightType::Direction){
                DirectionalLight* dl = dynamic_cast<DirectionalLight*>(l.get());
                ldir = (dl->dir).normalize(); 
                att = 1.0f;
            } else { blinn_phong = blinn_phong + l.get()->color.toVec3() * l.get()->intensity; continue; }

            // Check if there's geometry between the fragment and the light source
            if(!((tri.flags | (debug & 0xffu)) & DISABLE_SHADOWS)){
                Vec3 sori = Vec3::dot(ldir, hit.normal) < 0.0f ? hit.point - hit.normal*EPSILON : hit.point + hit.normal*EPSILON;
                Ray lray = Ray(sori, ldir);
                Hit hit2;
                if(intersection_shader(lray, hit2, DISABLE_SHADOWS | DISABLE_SHADING)) continue;
            }

            // Compute fragment's specular and diffuse components
            float diff = max(0.0f, Vec3::dot(hit.phong, ldir));
            Vec3 diffuse = (l.get()->color.toVec3() * l.get()->intensity) * (diff * mat.diff);

            half = (ldir + view).normalize();
            float spec = powf(max(0.0f, Vec3::dot(hit.phong, half)), mat.spec);
            Vec3 specular = (l.get()->color.toVec3() * l.get()->intensity) * spec;

            blinn_phong = blinn_phong + (diffuse + specular) * att;
        }

    } else blinn_phong = Vec3(1.0f,1.0f,1.0f);

    // Compute final fragment color
    out = Fragment(blinn_phong, 1.0f) * tex;
    return out;
}

// Texture shader: compute texel from tri texture and hit UV coordinates
Fragment PolyRenderer::texture_shader(Hit& hit){
    Material& mat = mats[tris[hit.tri].mat];
    if(mat.texture){
        uint16_t tx = hit.u * (TEXTURE_SIZE - 1); tx %= (TEXTURE_SIZE-1);
        uint16_t ty = hit.v * (TEXTURE_SIZE - 1); ty %= (TEXTURE_SIZE-1);
        return Fragment(mat.texture[tx + ty*TEXTURE_SIZE]);
    } else return Fragment(mat.color);
}

// TODO Bump shader: compute surface normal from tri normal map and hit UV coordinates
Vec3 PolyRenderer::bump_shader(Hit& hit){
    Material& mat = mats[tris[hit.tri].mat];
    if(mat.bump){
        // TODO: convert to tangent coordinates and compute normal vector using bump texture
        uint16_t tx = hit.u * (TEXTURE_SIZE - 1); tx %= (TEXTURE_SIZE-1);
        uint16_t ty = hit.v * (TEXTURE_SIZE - 1); ty %= (TEXTURE_SIZE-1);
        return mat.bump[tx + ty*TEXTURE_SIZE].toVec3() * 2.0f - 1.0f;
    } else return Vec3(1.0f,1.0f,1.0f);
}

// Reflection shader: compute N reflection
Fragment PolyRenderer::reflection_shader(Ray& ray, Hit& hit, uint8_t N_REFLECTION){
    Material& mat = mats[tris[hit.tri].mat];
    Fragment frag = fragment_shader(hit);
    if(N_REFLECTION<MAX_REFLECTIONS){   // Check if we hit the limit of reflection rays
        Vec3 rdir = ray.dir - hit.phong * (Vec3::dot(ray.dir, hit.phong) * 2.0f);
        Ray rray = Ray(hit.point, rdir);
        Hit rhit;
        return intersection_shader(rray, rhit) ? frag * (1.0f-mat.reflective) + reflection_shader(rray, rhit, N_REFLECTION+1) * mat.reflective : frag;
    } else return frag;
}

// TODO Refraction shader: compute N refraction
Fragment PolyRenderer::refraction_shader(Ray& ray, Hit& hit, uint8_t N_REFRACTION){
    return Fragment();
}

// Save scene into .png file
void PolyRenderer::save(const char* path){
    savePNG(path, frame);
    printf("\e[1;93m render saved in '\e[95m%s\e[93m'\e[0m\n", path);
}

// Get Cpu code
const char* PolyRenderer::getCpu(){
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

// Print poly-classic intro
void PolyRenderer::printIntro(){
    system("clear");
    printf(" \e[1;91m▄▄▄   \e[92m▄▄   \e[94m▄  \e[95m▄   ▄ \n");
    printf(" \e[91m█  █ \e[92m█  █  \e[94m█   \e[95m█ █  \e[93m ▄▄ ▄   ▄   ▄▄  ▄▄ ▄  ▄▄\n");
    printf(" \e[91m█▀▀  \e[92m█  █  \e[94m█    \e[95m█  \e[93m █   █  █▄█ ▀▄  ▀▄  █ █\n");
    printf(" \e[91m█    \e[92m▀▄▄▀  \e[94m█▄▄  \e[95m█   \e[93m▀▄▄ █▄ █ █ ▄▄▀ ▄▄▀ █ ▀▄▄\n");
    printf("\e[91m         - diegojose.parragan@um.es -\n\e[0m\n");
}

// Print message to screen and save to internal string vector
void PolyRenderer::PolyRenderer::polyMsg(string msg){
    printf("%s", msg.c_str());
    msgvector.push_back(msg);
}

/* Load .png into RGBA buffer */
RGBA* PolyRenderer::loadPNG(const char* path){
FILE* input = fopen(path, "rb");
        if(!input){ PolyRenderer::polyMsg("\e[1;91m  err loading '" + string(path) + "': file could not be opened!\n\e[0m"); exit(EXIT_FAILURE); }

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if(!png){ 
            fclose(input);
            PolyRenderer::polyMsg("\e[1;91m  err texture '" + string(path) + "' read struct could not be loaded!\n\e[0m"); exit(EXIT_FAILURE);
        }

    png_infop info = png_create_info_struct(png);
        if(!info){
            fclose(input);
            png_destroy_read_struct(&png, NULL, NULL);
            PolyRenderer::polyMsg("\e[1;91m  err texture '" + string(path) + "' info struct could not be loaded!\e[0m\n"); exit(EXIT_FAILURE);
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

    if(fclose(input)==-1){ PolyRenderer::polyMsg("\e[1;91m  err closing '" + string(path) + "'\e[0m\n"); exit(EXIT_FAILURE); }
    PolyRenderer::polyMsg("\e[1;93m  loaded texture '" + string(path) + "'\e[0m\n");
    
    return texture;
}

/* Save RGBA buffer as .png */
void PolyRenderer::savePNG(const char* path, RGBA* texture){
FILE* output = fopen(path, "wb");
        if(!output){ PolyRenderer::polyMsg("\e[1;91m  err saving to '" + string(path) + "': file could not be opened!\n\e[0m"); exit(EXIT_FAILURE); }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if(!png){ 
            fclose(output);
            PolyRenderer::polyMsg("\e[1;91m  err saving render\n\e[0m");
            exit(EXIT_FAILURE);
        }

    png_infop info = png_create_info_struct(png);
        if(!info){ 
            fclose(output);
            png_destroy_write_struct(&png, NULL);
            PolyRenderer::polyMsg("\e[1;91m  err saving render\e[0m\n");
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

    if(fclose(output)==-1){ PolyRenderer::polyMsg("\e[1;91m  err closing '" + string(path) + "'\e[0m\n"); exit(EXIT_FAILURE); }
}
