/*
    POLY classic ~ main 
    Diego Párraga Nicolás ~ diegojose.parragan@um.es
*/

#include <iostream>
#include <cstdlib>
#include "PolyRenderer.h"

// Global variables
PolyRenderer* renderer;
string output;

/* Print usage */
void printUsage(){
    printf("\e[1;93m Usage: \e[91mp\e[92mo\e[94ml\e[95my \e[93mclassic \e[95m-i=SCENE.POLY \e[96m [-t=THREADS] [-o=OUTPUT_PATH]\e[0m\n");
}

void sigintHandler(int sig){
    printf("\n\e[1;96m SIGNINT received, rendering stopped\e[0m\n");
    renderer->save(output.c_str());
    exit(EXIT_FAILURE);
}

/* Application entrypoint */
int main(int argc, char** argv){

    // Parse arguments
    if(argc<2 || argc>4){ printUsage(); return EXIT_FAILURE; }

    bool ENABLE_RENDERING_WINDOW = false;
    uint8_t nthreads = static_cast<uint8_t>(omp_get_max_threads());
    string input, outputname;
    for(uint8_t i=1; i<argc; i++){
        string arg = string(argv[i]);
        if(arg.length()>3 && arg[0]=='-' && arg[2]=='='){
            string value = arg.substr(3, arg.length()-3);
            switch (arg[1]){
                case 't' :  // Number of threads. Defaults to omp_max_threads()
                    nthreads = stoi(value);
                    break;
                case 'o' :  // Output path
                    output = value + "/";
                    break;  
                case 'i' :  // Input scene script
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
    PolyRenderer::printIntro();

    // Create rendering struct and start render
    renderer = new PolyRenderer();
    signal(SIGINT, sigintHandler);
    if(renderer->loadScene(input.c_str()) && renderer->render(nthreads, ENABLE_RENDERING_WINDOW))
        renderer->save(output.c_str());

    free(renderer);

    return EXIT_SUCCESS;
}