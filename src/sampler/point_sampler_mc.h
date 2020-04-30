#pragma once

#include <cmath>

#include "sampler/value_sampler.h"
#include "vec.h"

class PointSamplerLambertian
{

private:
    ValueSampler<double> sampler_ = ValueSampler<double>(0, 1);

public:
    PointSamplerLambertian() {}
    Vec sample() const
    {
        double r1 = sampler_.sample();
        double r2 = sampler_.sample();
        double phi = 2 * M_PI * r1;
        double x = std::cos(phi) * std::sqrt(r2);
        double y = std::sin(phi) * std::sqrt(r2);
        double z = std::sqrt(1 - r2);
        return Vec(x, y, z);
    }
};

class PointSamplerToSphere
{
private:
    ValueSampler<double> sampler_ = ValueSampler<double>(0, 1);

public:
    PointSamplerToSphere() {}
    Vec sample(double radius, double distance_squared) const
    {
        double r1 = sampler_.sample();
        double r2 = sampler_.sample();
        double z = 1 + r2 * (std::sqrt(1 - radius * radius / distance_squared) - 1);
        double phi = 2 * M_PI * r1;
        double r3 = std::sqrt(1 - z * z);
        double x = std::cos(phi) * r3;
        double y = std::sin(phi) * r3;
        return Vec(x, y, z);
    }
};
