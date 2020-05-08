#include <iostream>

#include "photonmap.h"
#include "photonmap_sample.h"
#include "ppm_sample.h"
#include "progressive_photonmapping.h"
#include "render.h"
#include "stochastic_ppm.h"

int main(int argc, char **argv)
{
    std::cout << "Path tracing renderer: edupt base raytracer" << std::endl << std::endl;

    // 640x480の画像、(2x2) * 4 sample / pixel
    // edupt::render(640, 480, 4, 2);
    int width = 640;
    int height = 480;
    int samples = 1;
    int supersamples = 1;
    int photon_num = 50000;
    int gather_photon_radius = 32;
    int gahter_max_photon_num = 64;
    std::string scene_name = "sceneforppm3";

    int photonmap_num = 10;
    std::stringstream ss;

// return    ppm_sample::render();

    //  photon_num = 1000000;
    //  ss << "image_scene6_" << width << "_" << height << "_" << samples << "_" << supersamples << "_" << photon_num <<
    //  "_"
    //     << gather_photon_radius << "_" << gahter_max_photon_num << ".hdr";
    //  return photonmap::standard::render(ss.str(), width, height, samples, supersamples, photon_num,
    //  gather_photon_radius,
    //                                     gahter_max_photon_num);

    // ss << "image_scene6_multiple_" << width << "_" << height << "_" << samples << "_" << supersamples << "_" <<
    // photon_num << "_"
    //    << gather_photon_radius << "_" << gahter_max_photon_num << "_" << photonmap_num << ".hdr";
    // return photonmap::standard::render_multiple_photonmap(ss.str(), width, height, samples, supersamples, photon_num,
    //                                                       gather_photon_radius, gahter_max_photon_num,
    //                                                       photonmap_num);

    //  return photonmap_sample::render();

    samples = 4;
    supersamples = 4;
    gather_photon_radius = 100000;
    photonmap::progressive::InitialRadius = 25;
    gahter_max_photon_num = 32;
    photon_num = 10000000;
    ss << "image_" << scene_name << "_progressive7_" << width << "_" << height << "_" << samples << "_" << supersamples
       << "_" << photon_num << "_" << gather_photon_radius << "_" << gahter_max_photon_num << "_" << photonmap_num;

    return photonmap::progressive::render(ss.str(), width, height, samples, supersamples, photon_num,
                                          gather_photon_radius, gahter_max_photon_num);

//   return photonmap::sppm::render(ss.str(), width, height, samples, supersamples, photon_num, gather_photon_radius,
 //                                  gahter_max_photon_num);
}