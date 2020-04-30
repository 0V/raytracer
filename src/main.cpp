#include <iostream>

#include "photonmap.h"
#include "photonmap_sample.h"
#include "render.h"


int main(int argc, char **argv)
{
    std::cout << "Path tracing renderer: edupt" << std::endl << std::endl;

    // 640x480の画像、(2x2) * 4 sample / pixel
    // edupt::render(640, 480, 4, 2);
    return photonmap::render(640, 480, 2, 2, 500000, 32, 128);
    //    return photonmap_sample::render();
}