/*
    POLY classic ~ main 
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include <iostream>
#include <cstdlib>
#include "poly-classic.h"
using namespace std;

/* Print usage */
void printUsage(){
    printf("\e[1;93m Usage: \e[91mp\e[92mo\e[94ml\e[95my \e[93mclassic \e[95mSCENE.POLY \e[96m [THREADS]\e[0m\n");
}


/* Main function */
int main(int argc, char** argv){

    // Check arguments
    if(argc<2 || argc>3){ printUsage(); return EXIT_FAILURE; }
    uint NTHREADS = argc>2 ? static_cast<uint>(atoi(argv[2])) : static_cast<uint>(omp_get_max_threads());

    printIntro();

    filesystem::path scene(argv[1]);
    if(!filesystem::exists(scene)){ printf("\e[1;91m err loading scene '%s': file does not exist!\e[0m\n", argv[1]); return EXIT_FAILURE; }

    string name = scene.stem(), output = name + ".png";

    // Render
    PolyRenderer renderer = PolyRenderer();
    if(renderer.loadScene(argv[1]) && renderer.render(NTHREADS))
        renderer.save(output.c_str());

    return EXIT_SUCCESS;
}