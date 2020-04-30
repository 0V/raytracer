#pragma once

#include <cmath>

#include <random>
#include "sampler/value_sampler.h"
#include "vec.h"

class PointSamplerSphere
{
private:
    ValueSampler<double> sampler_ = ValueSampler<double>(-1, 1);

public:
    PointSamplerSphere() {}
    PointSamplerSphere(const double &radius_) : radius(radius_), square_radius(radius_ * radius_) {}
    const double radius = 1;
    const double square_radius = 1;

    Vec sample() const
    {
        Vec p;
        do
        {
            p = Vec(radius * sampler_.sample(), radius * sampler_.sample(), radius * sampler_.sample());
        } while (p.length_squared() >= square_radius);
        return p;
    }
};

class PointSamplerUnitSphere
{
private:
    ValueSampler<double> sampler_ = ValueSampler<double>(-1, 1);

public:
    PointSamplerUnitSphere() {}
    Vec sample() const
    {
        Vec p;
        do
        {
            p = Vec(sampler_.sample(), sampler_.sample(), sampler_.sample());
        } while (p.length_squared() >= 1);
        return p;
    }
};
