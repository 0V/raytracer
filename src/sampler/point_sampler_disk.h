#pragma once

#include <cmath>
#include <random>

#include "sampler/value_sampler.h"
#include "vec.h"

class PointSamplerDisk
{

private:
    ValueSampler<double> sampler_ = ValueSampler<double>(-1, 1);

public:
    PointSamplerDisk() {}
    PointSamplerDisk(const double &radius_) : radius(radius_), square_radius(radius_ * radius_) {}
    const double radius = 1;
    const double square_radius = 1;

    Vec sample() const
    {
        Vec p;
        do
        {
            p = Vec(radius * sampler_.sample(), radius * sampler_.sample(), 0);
        } while (p.length_squared() >= square_radius);
        return p;
    }
};

class PointSamplerUnitDisk
{
private:
    ValueSampler<double> sampler_ = ValueSampler<double>(-1, 1);

public:
    PointSamplerUnitDisk() {}

    Vec sample() const
    {
        Vec p;
        do
        {
            p = Vec(sampler_.sample(), sampler_.sample(), 0);
        } while (p.length_squared() >= 1);
        return p;
    }
};