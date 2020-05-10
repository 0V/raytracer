#pragma once

#include "ray.h"
#include "sampler/point_sampler_disk.h"

namespace edupt
{
    // Model : Thin Lens Model
    class DoFCamera
    {
    private:
        PointSamplerDisk sampler_;

    public:
        const double image_width;
        const double image_height;

        // ワールド座標系でのスクリーンの大きさ
        const double screen_width;
        const double screen_height;
        // スクリーンまでの距離
        const double screen_dist = 40.0;

        const Vec camera_position;
        const Vec camera_dir;
        const Vec camera_up;

        double supersamples;

        // スクリーンを張るベクトル
        const Vec screen_x;
        const Vec screen_y;
        const Vec screen_center;

        DoFCamera(const double width, const double height, const double screen_height_, const double screen_dist_,
                  const Vec& camera_position_, const Vec& camera_dir_, const Vec& camera_up_,
                  const double supersamples_, const ValueSampler<double>& sampler)
            : image_width(width),
              image_height(height),
              screen_width(screen_height_ * width / height),
              screen_height(screen_height_),
              screen_dist(screen_dist_),
              camera_position(camera_position_),
              camera_dir(camera_dir_),
              camera_up(camera_up_),
              supersamples(supersamples_),
              screen_x(normalize(cross(camera_dir, camera_up)) * screen_width),
              screen_y(normalize(cross(screen_x, camera_dir)) * screen_height),
              screen_center(camera_position + camera_dir * screen_dist)
        {
        }

        Ray get_ray(const double x, const double y, const int supersample_x = 0, const int supersample_y = 0)
        {
            const double rate = (1.0 / supersamples);
            const double r1 = supersample_x * rate + rate / 2.0;
            const double r2 = supersample_y * rate + rate / 2.0;
            const Vec screen_position =
                screen_center + screen_x * ((r1 + x) / image_width - 0.5) + screen_y * ((r2 + y) / image_height - 0.5);

            // レイを飛ばす方向
            const Vec dir = normalize(screen_position - camera_position);

            return Ray(camera_position, dir);
        }
    };
}  // namespace edupt