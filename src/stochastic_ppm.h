#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

#include "hdr.h"
#include "intersection.h"
#include "kdtree.h"
#include "material.h"
#include "pinhole_camera.h"
#include "dof_camera.h"
#include "ppm.h"
#include "radiance.h"
#include "random.h"
#include "ray.h"
#include "sampler/value_sampler.h"
#include "scene.h"
#include "vec.h"

namespace photonmap
{
    namespace stochasticppm
    {
        using namespace edupt;
        using namespace photonmap::utility;

        // クラスにしてフィールドにしたい
        double InitialRadius = 25;
        double Alpha = 0.7;  // the alpha parameter of PPM

        struct ProgressiveIntersection
        {
            Color accumulated_flux;
            Color emission;

            Intersection intersection;
            Vec position;
            double weight;
            double photon_radius_2;  // photon_radius^2
            int photon_count;
            int index;
            ProgressiveIntersection(Intersection intersection_, double weight_, double photon_radius_2_,
                                    int photon_count_, int index_, const Color& emission_)
                : intersection(intersection_),
                  position(intersection_.hitpoint.position),
                  weight(weight_),
                  photon_radius_2(photon_radius_2_),
                  index(index_),
                  emission(emission_)
            {
            }
        };

        struct Photon
        {
            Vec position;
            Color power;
            Vec incident;
            Photon(const Vec& position_, const Color& power_, const Vec& incident_)
                : position(position_), power(power_), incident(incident_)
            {
            }
        };

        using PointMap = KDTree<ProgressiveIntersection>;

        const int LightID = 0;
        const double INF = 1e20;
        const double EPS = 1e-6;
        const double MaxDepth = 5;

        // コサイン項によるImportance Sampling
        decltype(auto) cosine_sampling(const Vec& w, const Vec& u, const Vec& v, ValueSampler<double>& sampler01)
        {
            const double u1 = 2 * M_PI * sampler01.sample();
            const double u2 = sampler01.sample();
            const double u3 = std::sqrt(u2);
            return normalize((u * std::cos(u1) * u3 + v * std::sin(u1) * u3 + w * std::sqrt(1.0 - u2)));
        }

        // 球の１点をサンプリング
        decltype(auto) sphere_sampling(ValueSampler<double>& sampler01)
        {
            const double r1 = 2 * M_PI * sampler01.sample();
            const double r2 = 1.0 - 2.0 * sampler01.sample();
            const double r3 = 1.0 - r2 * r2;
            return Vec(std::sqrt(r3) * std::cos(r1), std::sqrt(r3) * std::sin(r1), r2);
        }

        void create_photon(const Sphere& light_sphere, ValueSampler<double>& sampler01, PointMap* point_map,
                           double gather_radius, double gahter_max_photon_num)
        {
            const Vec source_pos = light_sphere.position + ((light_sphere.radius + EPS) * sphere_sampling(sampler01));
            const Vec source_dir = normalize(source_pos - light_sphere.position);

            // 光源上の点から半球サンプリングする
            Vec w, u, v;
            w = source_dir;

            if (fabs(w.x) > 0.1)
            {
                u = normalize(cross(Vec(0.0, 1.0, 0.0), w));
            }
            else
            {
                u = normalize(cross(Vec(1.0, 0.0, 0.0), w));
            }
            v = cross(w, u);

            Vec light_dir = cosine_sampling(w, u, v, sampler01);

            Ray current_ray(source_pos, light_dir);

            Color current_flux = light_sphere.emission * 4.0 * std::pow(light_sphere.radius * M_PI, 2.0);

            bool trace_end = false;

            while (!trace_end)
            {
                if (std::max(current_flux.x, std::max(current_flux.y, current_flux.z)) <= 0.0) break;

                Intersection intersect_data;
                if (!intersect_scene(current_ray, &intersect_data)) break;
                const Sphere& obj = spheres[intersect_data.object_id];
                // 物体に対する入射方向が表か裏かを確認する
                const Vec orienting_normal = dot(intersect_data.hitpoint.normal, current_ray.dir) < 0.0
                                                 ? intersect_data.hitpoint.normal            // 物体外から入射
                                                 : (-1.0 * intersect_data.hitpoint.normal);  // 物体中から入射

                switch (obj.reflection_type)
                {
                    case edupt::REFLECTION_TYPE_DIFFUSE:
                    {
                        PointMap::ResultQueue result_queue;

                        PointMap::Query query(intersect_data.hitpoint.position, orienting_normal, gather_radius,
                                              gahter_max_photon_num);
                        point_map->SearchNearest(&result_queue, query);
                        // キューからフォトンを取り出しvectorに格納する
                        std::vector<PointMap::ElementForQueue> points;
                        points.reserve(result_queue.size());

                        //                        if (result_queue.size() > 0) std::cout << result_queue.size();
                        while (!result_queue.empty())
                        {
                            PointMap::ElementForQueue p = result_queue.top();
                            result_queue.pop();
                            points.push_back(p);
                            //                            max_distance2 = std::max(max_distance2, p.distance2);
                        }

                        for (auto& point : points)
                        {
                            if ((dot(point.point->intersection.hitpoint.normal, intersect_data.hitpoint.normal) > 1e-3))
                            {
                                if (point.point->photon_count == 0)
                                {
                                    point.point->photon_radius_2 =
                                        (intersect_data.hitpoint.position - point.point->position).length_squared();
                                    point.point->accumulated_flux =
                                        point.point->accumulated_flux + (multiply(obj.color, current_flux) / M_PI);

                                    point.point->photon_count++;
                                }
                                else
                                {
                                    if (point.distance2 < point.point->photon_radius_2)
                                    {
                                        const int N = point.point->photon_count;
                                        const int M = 1;

                                        double g2 = (N + Alpha * M) / (N + M);
                                        point.point->photon_radius_2 = point.point->photon_radius_2 * g2;
                                        point.point->photon_count++;
                                        // point.point->accumulated_flux =
                                        //     obj.emission + (point.point->accumulated_flux +
                                        //                     multiply(point.point->weight, current_flux) / M_PI) *
                                        //                        g;
                                        point.point->accumulated_flux = (point.point->accumulated_flux +
                                                                         (multiply(obj.color, current_flux) / M_PI)) *
                                                                        g2;
                                    }
                                }
                            }
                        }

                        // A Practical Guide to Global Illumination using Photon Mapsとは異なるが
                        // RGBの平均値を反射確率とする。
                        // TODO: Depthに応じて上げたい
                        const double probability = (obj.color.x + obj.color.y + obj.color.z) / 3;
                        if (probability > sampler01.sample())
                        {
                            // 反射
                            // orienting_normalの方向を基準とした正規直交基底(w, u,
                            // v)を作り。この基底に対する半球内で次のレイを飛ばす。
                            Vec diffuse_w, diffuse_u, diffuse_v;
                            diffuse_w = orienting_normal;
                            if (fabs(diffuse_w.x) > 0.1)
                            {
                                diffuse_u = normalize(cross(Vec(0.0, 1.0, 0.0), diffuse_w));
                            }
                            else
                            {
                                diffuse_u = normalize(cross(Vec(1.0, 0.0, 0.0), diffuse_w));
                            }
                            diffuse_v = cross(diffuse_w, diffuse_u);
                            Vec dir = cosine_sampling(diffuse_w, diffuse_u, diffuse_v, sampler01);

                            current_ray = Ray(intersect_data.hitpoint.position, dir);
                            current_flux = multiply(current_flux, obj.color) / probability;
                        }
                        else
                        {  // 吸収
                            trace_end = true;
                        }
                    }
                    break;
                    case edupt::REFLECTION_TYPE_SPECULAR:
                    {
                        // 完全鏡面
                        current_ray = Ray(intersect_data.hitpoint.position,
                                          reflection(current_ray.dir, intersect_data.hitpoint.normal));
                        current_flux = multiply(current_flux, obj.color);
                    }
                    break;
                    case edupt::REFLECTION_TYPE_REFRACTION:
                    {
                        // 屈折
                        Ray reflection_ray = Ray(intersect_data.hitpoint.position,
                                                 reflection(current_ray.dir, intersect_data.hitpoint.normal));
                        // レイの屈折方向がオブジェクトの内側方向か外側方向かを確認する
                        const bool is_into = dot(intersect_data.hitpoint.normal, orienting_normal) >
                                             0.0;  // レイがオブジェクトから出るのか、入るのか

                        // Snellの法則
                        const double nc = 1.0;  // 真空の屈折率
                        // edupt::kIor オブジェクトの屈折率
                        const double nt = edupt::kIor;  // オブジェクトの屈折率
                        const double nnt = is_into ? nc / nt : nt / nc;
                        const double ddn = dot(current_ray.dir, orienting_normal);
                        const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

                        if (cos2t < 0.0)
                        {  // 全反射
                            current_ray = reflection_ray;
                            current_flux = multiply(current_flux, obj.color);
                            continue;
                        }
                        // 屈折していく方向
                        Vec tdir =
                            normalize(current_ray.dir * nnt - intersect_data.hitpoint.normal * (is_into ? 1.0 : -1.0) *
                                                                  (ddn * nnt + std::sqrt(cos2t)));

                        // SchlickによるFresnelの反射係数の近似
                        const double a = nt - nc, b = nt + nc;
                        const double R0 = (a * a) / (b * b);
                        const double c = 1.0 - (is_into ? -ddn : dot(tdir, intersect_data.hitpoint.normal));
                        const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
                        const double Tr = 1.0 - Re;  // 屈折光の運ぶ光の量
                        const double probability = Re;

                        // 屈折と反射のどちらか一方を追跡する。
                        // ロシアンルーレットで決定する。
                        if (sampler01.sample() < probability)
                        {  // 反射
                            current_ray = reflection_ray;
                            // Fresnel係数Reを乗算し、ロシアンルーレット確率prob.で割る。
                            // 今、prob.=Reなので Re / prob. = 1.0 となる。
                            // よって、current_flux = Multiply(current_flux, obj.color) * Re / probability;
                            // が以下の式になる。 屈折の場合も同様。
                            current_flux = multiply(current_flux, obj.color);
                            continue;
                        }
                        else
                        {  // 屈折
                            current_ray = Ray(intersect_data.hitpoint.position, tdir);
                            current_flux = multiply(current_flux, obj.color);
                            continue;
                        }
                    }

                    break;
                }
            }
        }

        void create_point(const Ray& ray, ValueSampler<double>& sampler01, const int depth, double weight,
                          PointMap* point_map, int index)
        {
            Intersection intersection;
            // シーンと交差判定
            if (!intersect_scene(ray, &intersection)) return;

            const Sphere& now_object = spheres[intersection.object_id];
            const Hitpoint& hitpoint = intersection.hitpoint;
            const Vec orienting_normal =
                dot(hitpoint.normal, ray.dir) < 0.0 ? hitpoint.normal : (-1.0 * hitpoint.normal);

            // 交差位置の法線（物体からのレイの入出を考慮）
            // 色の反射率最大のものを得る。ロシアンルーレットで使う。
            // ロシアンルーレットの閾値は任意だが色の反射率等を使うとより良い。
            double russian_roulette_probability =
                std::max(now_object.color.x, std::max(now_object.color.y, now_object.color.z));

            // 反射回数が一定以上になったらロシアンルーレットの確率を急上昇させる。（スタックオーバーフロー対策）
            if (depth > kDepthLimit) russian_roulette_probability *= std::pow(0.5, depth - kDepthLimit);

            // ロシアンルーレットを実行し追跡を打ち切るかどうかを判断する。
            // ただしDepth回の追跡は保障する。
            if (depth > kDepth)
            {
                if (sampler01.sample() >= russian_roulette_probability) return;
            }
            else
                russian_roulette_probability = 1.0;  // ロシアンルーレット実行しなかった

            switch (now_object.reflection_type)
            {
                // 完全拡散面
                case REFLECTION_TYPE_DIFFUSE:
                {
                    point_map->AddData(ProgressiveIntersection(intersection, weight / russian_roulette_probability,
                                                               InitialRadius, 0, index, now_object.emission));
                    return;
                }
                break;
                // 完全鏡面
                case REFLECTION_TYPE_SPECULAR:
                {
                    // 完全鏡面なのでレイの反射方向は決定的。
                    // ロシアンルーレットの確率で除算するのは上と同じ。
                    create_point(Ray(hitpoint.position, reflection(ray.dir, hitpoint.normal)), sampler01, depth + 1,
                                 weight / russian_roulette_probability, point_map, index);
                    return;
                }
                break;

                // 屈折率kIorのガラス
                case REFLECTION_TYPE_REFRACTION:
                {
                    const Ray reflection_ray = Ray(hitpoint.position, reflection(ray.dir, hitpoint.normal));
                    const bool into =
                        dot(hitpoint.normal, orienting_normal) > 0.0;  // レイがオブジェクトから出るのか、入るのか

                    // Snellの法則
                    const double nc = 1.0;   // 真空の屈折率
                    const double nt = kIor;  // オブジェクトの屈折率
                    const double nnt = into ? nc / nt : nt / nc;
                    const double ddn = dot(ray.dir, orienting_normal);
                    const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

                    if (cos2t < 0.0)
                    {  // 全反射
                        create_point(reflection_ray, sampler01, depth + 1, weight / russian_roulette_probability,
                                     point_map, index);
                        return;
                    }
                    // 屈折の方向
                    const Ray refraction_ray =
                        Ray(hitpoint.position, normalize(ray.dir * nnt - hitpoint.normal * (into ? 1.0 : -1.0) *
                                                                             (ddn * nnt + std::sqrt(cos2t))));

                    // SchlickによるFresnelの反射係数の近似
                    const double a = nt - nc, b = nt + nc;
                    const double R0 = (a * a) / (b * b);

                    const double c = 1.0 - (into ? -ddn : dot(refraction_ray.dir, hitpoint.normal));
                    // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
                    const double Re = R0 + (1.0 - R0) * std::pow(c, 5.0);
                    // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
                    const double nnt2 = std::pow(into ? nc / nt : nt / nc, 2.0);
                    const double Tr = (1.0 - Re) * nnt2;  // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

                    // 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
                    // ロシアンルーレットで決定する。
                    const double probability = 0.25 + 0.5 * Re;
                    if (depth > 2)
                    {
                        if (sampler01.sample() < probability)
                        {  // 反射

                            create_point(reflection_ray, sampler01, depth + 1,
                                         weight * Re / (probability * russian_roulette_probability), point_map, index);
                            return;
                        }
                        else
                        {  // 屈折

                            create_point(refraction_ray, sampler01, depth + 1,
                                         weight * Tr / ((1.0 - probability) * russian_roulette_probability), point_map,
                                         index);
                            return;
                        }
                    }
                    else
                    {  // 屈折と反射の両方を追跡
                        create_point(reflection_ray, sampler01, depth + 1, weight * Re / russian_roulette_probability,
                                     point_map, index);
                        create_point(refraction_ray, sampler01, depth + 1, weight * Tr / russian_roulette_probability,
                                     point_map, index);
                        return;
                    }
                }
                break;
            }
        }

        void update_photons(int photon_count, PointMap* point_map, double gather_radius, double gahter_max_photon_num)
        {
            ValueSampler<double> sampler01(0, 1);
            //            double max_radius = 0;
            int devide_num = 100;
            int tmp_devide_num = photon_count;
            while (tmp_devide_num > 10000)
            {
                devide_num *= 100;
                tmp_devide_num /= 100;
                /* code */
            }

            for (size_t i = 0; i < photon_count; i++)
            {
                if (i % devide_num == 0)
                {
                    // double max_radius_2 = 0;
                    // for (const auto& p : point_map->GetData())
                    // {
                    //     max_radius_2 = std::max(max_radius_2, p.photon_radius_2);
                    // }
                    // max_radius = std::sqrt(max_radius_2);
                    std::cout << "Photon Tracing (i = " << i << ") " << (100.0 * i / (photon_count - 1)) << "%"
                              << std::endl;
                    //                    std::cout << "Rmax = " << max_radius << std::endl;
                }

                create_photon(spheres[LightID], sampler01, point_map, gather_radius, gahter_max_photon_num);
            }
        };

        void update_photons_with_render(int photon_count, PointMap* point_map, double gather_radius,
                                        double gahter_max_photon_num, int width, int height, double spp,
                                        std::string filename)
        {
            ValueSampler<double> sampler01(0, 1);
            double max_radius = 0;
            int devide_num = 1;
            int tmp_devide_num = photon_count;
            while (tmp_devide_num > 10000)
            {
                devide_num *= 100;
                tmp_devide_num /= 100;
                /* code */
            }

            Color* image = new Color[width * height];
            for (size_t i = 0; i < photon_count; i++)
            {
                create_photon(spheres[LightID], sampler01, point_map, gather_radius, gahter_max_photon_num);

                if (i % devide_num == 0)
                {
                    // double max_radius_2 = 0;
                    // for (const auto& p : point_map->GetData())
                    // {
                    //     max_radius_2 = std::max(max_radius_2, p.photon_radius_2);
                    // }
                    // max_radius = std::sqrt(max_radius_2);
                    std::cout << "Photon Tracing (i = " << i << ") " << (100.0 * i / (photon_count - 1)) << "%"
                              << std::endl;
                    //                    std::cout << "Rmax = " << max_radius << std::endl;

                    if (i % (devide_num * 10) == 0)
                    {
                        std::cout << "Rendering..." << std::endl;
                        for (size_t i = 0; i < width * height; i++)
                        {
                            image[i] = Color();
                        }
                        for (auto node : point_map->GetData())
                        {
                            // if (node.photon_count > 0)
                            // {
                            image[node.index] = image[node.index] +
                                                (node.emission + node.weight * node.accumulated_flux *
                                                                     (1.0 / (M_PI * node.photon_radius_2 * (i + 1)))) /
                                                    spp;
                            //                            }
                        }

                        save_ppm_file(filename + "_" + std::to_string(i) + ".ppm", image, width, height);
                        save_hdr_file(filename + "_" + std::to_string(i) + ".hdr", image, width, height);
                        //                    save_hdr_file(filename + "_" + std::to_string(i) + ".hdr", image, width,
                        //                    height);
                        std::cout << "Done Rendering" << std::endl;
                    }
                }
            }
        };

        int render(const std::string& filename, const int width, const int height, const int samples,
                   const int supersamples, int photon_num, double gather_photon_radius, int gahter_max_photon_num)
        {
            // カメラ位置
            const Vec camera_position = Vec(50.0, 52.0, 220.0);
            const Vec camera_dir = normalize(Vec(0.0, -0.04, -1.0));
            const Vec camera_up = Vec(0.0, 1.0, 0.0);

            // ワールド座標系でのスクリーンの大きさ
            const double screen_width = 30.0 * width / height;
            const double screen_height = 30.0;
            // スクリーンまでの距離
            const double screen_dist = 40.0;

            PinholeCamera camera(width, height, screen_height, screen_dist, camera_position, camera_dir, camera_up,
                                 supersamples);

            Color* image = new Color[width * height];
            std::vector<ProgressiveIntersection> hitpoint_list;
            int spp = samples * (supersamples * supersamples);

            std::cout << width << "x" << height << " " << spp << " spp" << std::endl;

            PointMap point_map;

            //            create_pointmap(width, height, &point_map, samples);

            // ValueSampler<double> sampler01(0, 1);
            // for (int y = 0; y < height; y++)
            // {
            //     std::cout << "Creating Point Tracing Map (y = " << y << ") " << (100.0 * y / (height - 1)) << "%"
            //               << std::endl;

            //     for (int x = 0; x < width; x++)
            //     {
            //         const int image_index = (height - y - 1) * width + x;
            //         // 一つのサブピクセルあたりsamples回サンプリングする
            //         for (int s = 0; s < samples; s++)
            //         {
            //             // スクリーン上の位置
            //             const Vec screen_position = screen_center + screen_x * ((double)x / width - 0.5) +
            //                                         screen_y * ((double)y / height - 0.5);
            //             // レイを飛ばす方向
            //             const Vec dir = normalize(screen_position - camera_position);

            //             create_point(Ray(camera_position, dir), sampler01, 0, 1, &point_map, image_index);
            //         }
            //     }
            // }

            ValueSampler<double> sampler01(0, 1);
            for (int y = 0; y < height; y++)
            {
                std::cout << "Creating Point Tracing Map (y = " << y << ") " << (100.0 * y / (height - 1)) << "%"
                          << std::endl;

                for (int x = 0; x < width; x++)
                {
                    const int image_index = (height - y - 1) * width + x;
                    for (int sy = 0; sy < supersamples; sy++)
                    {
                        for (int sx = 0; sx < supersamples; sx++)
                        {
                            // 一つのサブピクセルあたりsamples回サンプリングする
                            for (int s = 0; s < samples; s++)
                            {
                                create_point(camera.get_ray(x, y, sx, sy), sampler01, 0, 1.0, &point_map, image_index);
                            }
                        }
                    }
                }
            }

            std::cout << "Creating KDTree..." << std::endl;
            point_map.CreateKDtree();
            std::cout << "Done Creating KDTree" << std::endl;

            std::cout << "Emitting Photon..." << std::endl;

            update_photons_with_render(photon_num, &point_map, gather_photon_radius, gahter_max_photon_num, width,
                                       height, spp, filename);
            //            update_photons(photon_num, &point_map, gather_photon_radius, gahter_max_photon_num);
            std::cout << "Done Emitting Photon" << std::endl;

            std::cout << "Rendering..." << std::endl;
            for (auto node : point_map.GetData())
            {
                image[node.index] =
                    image[node.index] + (node.emission + node.weight * node.accumulated_flux *
                                                             (1.0 / (M_PI * node.photon_radius_2 * photon_num))) /
                                            spp;
            }
            std::cout << "Done Rendering" << std::endl;

            // 出力
            //        save_ppm_file(std::string("image.ppm"), image, width, height);
            save_ppm_file(filename + ".ppm", image, width, height);
            save_hdr_file(filename + ".hdr", image, width, height);
            return 0;
        }
    }  // namespace stochasticppm
}  // namespace photonmap