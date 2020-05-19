#include <iostream>

#include "photonmap.h"
#include "photonmap_sample.h"
#include "ppm_sample.h"
#include "progressive_photonmapping.h"
#include "raw_ppm.h"
#include "render.h"
#include "smallppm_exp.h"
#include "stochastic_ppm.h"

int main(int argc, char **argv)
{
    std::cout << "Path tracing renderer: edupt base raytracer" << std::endl << std::endl;

    // 640x480の画像、(2x2) * 4 sample / pixel
    // edupt::render(640, 480, 4, 2);
    constexpr int width = 320;
    constexpr int height = 240;
    int samples = 100;
    int supersamples = 4;
    int photon_num = 50000;
    int gather_photon_radius = 32;
    int gahter_max_photon_num = 64;
    std::string scene_name = "sceneforppm3";

    int photonmap_num = 10;
    std::stringstream ss;
 //   return edupt::render_dof(width, height, samples, supersamples, 180, 4);
    return edupt::render(width, height, samples, supersamples);

    // for (size_t i = 60; i < 200; i = i + 10)
    // {
    //    edupt::render_dof(width, height, samples, supersamples, i);
    // }
    // return 0;

    //    return photonmap::smallppm::render();

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

    //////////////// COMPACT SET ////////////////
    // samples = 32;
    // supersamples = 1;
    // gather_photon_radius = 25;
    // photonmap::progressive::InitialRadius = 25;
    // gahter_max_photon_num = 10000;
    // photon_num = 10000000;
    /////////////////////////////////////////////

    samples = 1;
    supersamples = 4;
    gather_photon_radius = 5;
    photonmap::progressive::InitialRadius = gather_photon_radius;

    gahter_max_photon_num = 1000000;
    photon_num = 100000;

    ss << "image_" << scene_name << "_ppm101_" << width << "_" << height << "_" << samples << "_" << supersamples << "_"
       << photon_num << "_" << gather_photon_radius << "_" << gahter_max_photon_num << "_" << photonmap_num;
   //  return photonmap::progressive2::render(ss.str(), width, height, samples, supersamples, photon_num,
   //                                         gather_photon_radius, gahter_max_photon_num);

   //  return photonmap::progressive::render(ss.str(), width, height, samples, supersamples, photon_num,
   //                                        gather_photon_radius, gahter_max_photon_num);

    photonmap::sppm::StochasticPpm<width, height> sppm(samples, supersamples);
    return sppm.render(ss.str(), width, height, samples, supersamples, photon_num, gather_photon_radius,
                       gahter_max_photon_num, 100);

    //   return photonmap::sppm::render(ss.str(), width, height, samples, supersamples, photon_num,
    //   gather_photon_radius,
    //                                  gahter_max_photon_num);
}