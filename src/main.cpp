#include <iostream>

#include "photonmap.h"
#include "photonmap_sample.h"
#include "progressive_photonmapping.h"
#include "render.h"

int main(int argc, char **argv)
{
    std::cout << "Path tracing renderer: edupt" << std::endl << std::endl;

    // 640x480の画像、(2x2) * 4 sample / pixel
    // edupt::render(640, 480, 4, 2);
    int width = 640;
    int height = 480;
    int samples = 2;
    int supersamples = 2;
    int photon_num = 50000;
    int gather_photon_radius = 32;
    int gahter_max_photon_num = 64;
    std::string scene_name = "sceneforppm2";

    int photonmap_num = 10;
    std::stringstream ss;


    ss << "image_scene6_" << width << "_" << height << "_" << samples << "_" << supersamples << "_" << photon_num <<
    "_"
       << gather_photon_radius << "_" << gahter_max_photon_num << ".hdr";
    return photonmap::standard::render(ss.str(), width, height, samples, supersamples, photon_num,
    gather_photon_radius,
                                       gahter_max_photon_num);

    // ss << "image_scene6_multiple_" << width << "_" << height << "_" << samples << "_" << supersamples << "_" <<
    // photon_num << "_"
    //    << gather_photon_radius << "_" << gahter_max_photon_num << "_" << photonmap_num << ".hdr";
    // return photonmap::standard::render_multiple_photonmap(ss.str(), width, height, samples, supersamples, photon_num,
    //                                                       gather_photon_radius, gahter_max_photon_num,
    //                                                       photonmap_num);

    // return photonmap_sample::render();

    samples = 1;
    supersamples = 2;
    gather_photon_radius = 200;
    photonmap::progressive::InitialRadius = gather_photon_radius;
    gahter_max_photon_num = 2000000;
    photon_num = 10000000;

    ss << "image_" << scene_name << "_progressive_" << width << "_" << height << "_" << samples << "_" << supersamples
       << "_" << photon_num << "_" << gather_photon_radius << "_" << gahter_max_photon_num << "_" << photonmap_num
       << ".hdr";
    return photonmap::progressive::render(ss.str(), width, height, samples, supersamples, photon_num,
                                          gather_photon_radius, gahter_max_photon_num);
}