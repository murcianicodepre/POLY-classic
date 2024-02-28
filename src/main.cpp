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
    printf("\e[1;93m Usage: \e[91mp\e[92mo\e[94ml\e[95my \e[93mclassic \e[95m-i=SCENE.POLY \e[96m [-t=THREADS] [-o=OUTPUT_PATH]\e[0m\n");
}

/* Application entrypoint */
int main(int argc, char** argv){

    // Parse arguments
    if(argc<2 || argc>4){ printUsage(); return EXIT_FAILURE; }

    uint nthreads = static_cast<uint8_t>(omp_get_max_threads());
    string input, output, outputname;
    for(uint i=1; i<argc; i++){
        string arg = string(argv[i]);
        if(arg.length()>3 && arg[0]=='-' && arg[2]=='='){
            string value = arg.substr(3, arg.length()-3);
            switch (arg[1]){
                case 't' :
                    nthreads = stoi(value);
                    break;
                case 'o' :
                    output = value + "/";
                    break;
                case 'i' :
                    input = value;
                    filesystem::path scene(value);
                    if(!filesystem::exists(scene)){ printf("\e[1;91m err loading scene '%s': file does not exist!\e[0m\n", value.c_str()); return EXIT_FAILURE; }
                    else if(scene.extension().string() != ".poly"){ printf("\e[1;91m err loading scene '%s': not a .poly file!\e[0m\n", value.c_str()); return EXIT_FAILURE; }
                    outputname = scene.stem(); outputname += ".png";
                    break;
            }
        } else { printUsage(); return EXIT_FAILURE; }
    }

    // Compute final output path
    output += outputname;

    // Print intro
    printIntro();

    // Create rendering struct and start render
    PolyRenderer renderer = PolyRenderer();
    if(renderer.loadScene(input.c_str()) && renderer.render(nthreads))
        renderer.save(output.c_str());

    return EXIT_SUCCESS;
}