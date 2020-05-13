#pragma once

#include <cmath>
#include <random>

#include "../vec.h"
#include "value_sampler.h"

class PointSamplerDisk
{
private:
    std::mt19937& engine_;
    std::uniform_real_distribution<double> dist_ = std::uniform_real_distribution<double>(-1, 1);

public:
    PointSamplerDisk() = default;
    PointSamplerDisk(std::mt19937& engine) : engine_(engine) {}
    PointSamplerDisk(std::mt19937& engine, const double& radius_)
        : engine_(engine), radius(radius_), square_radius(radius_ * radius_)
    {
    }

    const double radius = 1;
    const double square_radius = 1;

    void set_engine(std::mt19937& engine) { engine_ = engine; }

    Vec sample()
    {
        double r, theta;
        const double sx = dist_(engine_);
        const double sy = dist_(engine_);

        // sx, syが0,0だった場合は特別に処理
        if (sx == 0.0 && sy == 0.0)
        {
            return Vec(0, 0, 0);
        }
        // 四つに分割した円の各部位で別々の処理になる
        if (sx >= -sy)
        {
            if (sx > sy)
            {
                r = sx;
                if (sy > 0.0)
                    theta = sy / r;
                else
                    theta = 8.0 + sy / r;
            }
            else
            {
                r = sy;
                theta = 2.0 - sx / r;
            }
        }
        else
        {
            if (sx <= sy)
            {
                r = -sx;
                theta = 4.0 - sy / r;
            }
            else
            {
                r = -sy;
                theta = 6.0 + sx / r;
            }
        }
        theta *= M_PI / 4.0;
        return Vec(r * cosf(theta), r * sinf(theta), 0) * radius;
    }
};

class PointSamplerDisk_Const
{
private:
    mutable std::random_device seed_gen_;
    mutable std::mt19937 engine_ = std::mt19937(seed_gen_());
    mutable std::uniform_real_distribution<double> dist_ = std::uniform_real_distribution<double>(-1, 1);

public:
    PointSamplerDisk_Const() {}
    PointSamplerDisk_Const(const double& radius_) : radius(radius_), square_radius(radius_ * radius_) {}

    const double radius = 1;
    const double square_radius = 1;

    Vec sample() const
    {
        double r, theta;
        const double sx = dist_(engine_);
        const double sy = dist_(engine_);

        // sx, syが0,0だった場合は特別に処理
        if (sx == 0.0 && sy == 0.0)
        {
            return Vec(0, 0, 0);
        }
        // 四つに分割した円の各部位で別々の処理になる
        if (sx >= -sy)
        {
            if (sx > sy)
            {
                r = sx;
                if (sy > 0.0)
                    theta = sy / r;
                else
                    theta = 8.0 + sy / r;
            }
            else
            {
                r = sy;
                theta = 2.0 - sx / r;
            }
        }
        else
        {
            if (sx <= sy)
            {
                r = -sx;
                theta = 4.0 - sy / r;
            }
            else
            {
                r = -sy;
                theta = 6.0 + sx / r;
            }
        }
        theta *= M_PI / 4.0;
        return Vec(r * cosf(theta), r * sinf(theta), 0) * radius;
    }
};
