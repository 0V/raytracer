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
#include "kdtree.h"
#include "material.h"
#include "ppm.h"
#include "radiance.h"
#include "random.h"
#include "sampler/value_sampler.h"
#include "vec.h"

namespace photonmap
{
    using namespace edupt;

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

    using PhotonMap = KDTree<Photon>;

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

    // ray方向からの放射輝度を求める
    Color photonmap_radiance(const Ray& ray, ValueSampler<double>& sampler01, const int depth, PhotonMap* photon_map,
                             const double gather_radius, const int gahter_max_photon_num)
    {
        Intersection intersection;
        // シーンと交差判定
        if (!intersect_scene(ray, &intersection)) return kBackgroundColor;

        const Sphere& now_object = spheres[intersection.object_id];
        const Hitpoint& hitpoint = intersection.hitpoint;
        const Vec orienting_normal = dot(hitpoint.normal, ray.dir) < 0.0
                                         ? hitpoint.normal
                                         : (-1.0 * hitpoint.normal);  // 交差位置の法線（物体からのレイの入出を考慮）
        // 色の反射率最大のものを得る。ロシアンルーレットで使う。
        // ロシアンルーレットの閾値は任意だが色の反射率等を使うとより良い。
        double russian_roulette_probability =
            std::max(now_object.color.x, std::max(now_object.color.y, now_object.color.z));

        // 反射回数が一定以上になったらロシアンルーレットの確率を急上昇させる。（スタックオーバーフロー対策）
        if (depth > kDepthLimit) russian_roulette_probability *= pow(0.5, depth - kDepthLimit);

        // ロシアンルーレットを実行し追跡を打ち切るかどうかを判断する。
        // ただしDepth回の追跡は保障する。
        if (depth > kDepth)
        {
            if (sampler01.sample() >= russian_roulette_probability) return now_object.emission;
        }
        else
            russian_roulette_probability = 1.0;  // ロシアンルーレット実行しなかった

        Color incoming_radiance;

        switch (now_object.reflection_type)
        {
            // 完全拡散面
            case REFLECTION_TYPE_DIFFUSE:
            {
                PhotonMap::ResultQueue result_queue;

                PhotonMap::Query query(hitpoint.position, orienting_normal, gather_radius, gahter_max_photon_num);
                photon_map->SearchNearest(&result_queue, query);
                Color accumulated_flux;
                double max_distance2 = -1;

                // キューからフォトンを取り出しvectorに格納する
                std::vector<PhotonMap::ElementForQueue> photons;
                photons.reserve(result_queue.size());
                while (!result_queue.empty())
                {
                    PhotonMap::ElementForQueue p = result_queue.top();
                    result_queue.pop();
                    photons.push_back(p);
                    max_distance2 = std::max(max_distance2, p.distance2);
                }

                // Coneフィルタ
                const double max_distance = std::sqrt(max_distance2);
                const double k = 1.1;
                for (int i = 0; i < photons.size(); i++)
                {
                    const double w = 1.0 - (std::sqrt(photons[i].distance2) / (k * max_distance));  // 重み
                    const Color v =
                        multiply(now_object.color, photons[i].point->power) / M_PI;  // DiffuseのBRDF = 1.0 / π
                    accumulated_flux = accumulated_flux + w * v;
                }

                accumulated_flux = accumulated_flux / (1.0 - 2.0 / (3.0 * k));

                // // Gaussianフィルタ
                // const double max_distance = std::sqrt(max_distance2);
                // const double alpha = 0.918;
                // const double beta = 1.953;
                // for (int i = 0; i < photons.size(); i++)
                // {
                //     const double w =
                //         alpha *
                //         (1.0 - ((1 - std::exp(-beta * photons[i].distance2 * 0.5 / (max_distance * max_distance)) /
                //                          (1 - std::exp(-beta)))));  // 重み
                //     const Color v =
                //         multiply(now_object.color, photons[i].point->power) / M_PI;  // DiffuseのBRDF = 1.0 / π
                //     accumulated_flux = accumulated_flux + w * v;
                // }

                if (max_distance2 > 0.0)
                {
                    return now_object.emission +
                           accumulated_flux / (M_PI * max_distance2) / russian_roulette_probability;
                }
            }
            break;

            // 完全鏡面
            case REFLECTION_TYPE_SPECULAR:
            {
                // 完全鏡面なのでレイの反射方向は決定的。
                // ロシアンルーレットの確率で除算するのは上と同じ。
                return now_object.emission +
                       photonmap_radiance(
                           Ray(hitpoint.position, ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir)),
                           sampler01, depth + 1, photon_map, gather_radius, gahter_max_photon_num) /
                           russian_roulette_probability;
            }
            break;

            // 屈折率kIorのガラス
            case REFLECTION_TYPE_REFRACTION:
            {
                const Ray reflection_ray =
                    Ray(hitpoint.position, ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir));
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
                    return now_object.emission +
                           multiply(photonmap_radiance(reflection_ray, sampler01, depth + 1, photon_map, gather_radius,
                                                       gahter_max_photon_num),
                                    now_object.color) /
                               russian_roulette_probability;
                }
                // 屈折の方向
                const Ray refraction_ray =
                    Ray(hitpoint.position,
                        normalize(ray.dir * nnt - hitpoint.normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t))));

                // SchlickによるFresnelの反射係数の近似
                const double a = nt - nc, b = nt + nc;
                const double R0 = (a * a) / (b * b);

                const double c = 1.0 - (into ? -ddn : dot(refraction_ray.dir, -1.0 * orienting_normal));
                // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
                const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
                // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
                const double nnt2 = pow(into ? nc / nt : nt / nc, 2.0);
                const double Tr = (1.0 - Re) * nnt2;  // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

                // 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
                // ロシアンルーレットで決定する。
                const double probability = 0.25 + 0.5 * Re;
                if (depth > 2)
                {
                    if (sampler01.sample() < probability)
                    {  // 反射
                        return now_object.emission +
                               multiply(now_object.color,
                                        photonmap_radiance(reflection_ray, sampler01, depth + 1, photon_map,
                                                           gather_radius, gahter_max_photon_num) *
                                            Re) /
                                   (probability * russian_roulette_probability);
                    }
                    else
                    {  // 屈折
                        return now_object.emission +
                               multiply(now_object.color,
                                        photonmap_radiance(refraction_ray, sampler01, depth + 1, photon_map,
                                                           gather_radius, gahter_max_photon_num) *
                                            Re) /
                                   ((1.0 - probability) * russian_roulette_probability);
                    }
                }
                else
                {  // 屈折と反射の両方を追跡
                    return now_object.emission +
                           multiply(now_object.color,
                                    photonmap_radiance(reflection_ray, sampler01, depth + 1,
                                     photon_map, gather_radius,
                                                       gahter_max_photon_num) *
                                            Re +
                                        photonmap_radiance(refraction_ray, sampler01, depth + 1, photon_map,
                                                           gather_radius, gahter_max_photon_num) *
                                            Tr) /
                               russian_roulette_probability;
                    ;
                }
            }
            break;
        }

        return Color();
    }

    void create_photonmap(const int photon_count, PhotonMap* photon_map, const Sphere& light_sphere)
    {
        ValueSampler<double> sampler01(0, 1);
        for (int i = 0; i < photon_count; i++)
        {
            // 光源からフォトンを発射する
            // 光源は球。球の一点をサンプリングする
            const double r1 = 2 * M_PI * sampler01.sample();
            const double r2 = 1.0 - 2.0 * sampler01.sample();
            const double r3 = 1.0 - r2 * r2;

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

            Color current_flux = light_sphere.emission * 4.0 * std::pow(light_sphere.radius * M_PI, 2.0) / photon_count;

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
                        // 拡散面なのでフォトンをフォトンマップに格納する
                        photon_map->AddData(Photon(intersect_data.hitpoint.position, current_flux, current_ray.dir));

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
                                          current_ray.dir - intersect_data.hitpoint.normal * 2.0 *
                                                                dot(intersect_data.hitpoint.normal, current_ray.dir));
                        current_flux = multiply(current_flux, obj.color);
                    }
                    break;
                    case edupt::REFLECTION_TYPE_REFRACTION:
                    {
                        // 屈折
                        Ray reflection_ray =
                            Ray(intersect_data.hitpoint.position,
                                current_ray.dir - intersect_data.hitpoint.normal * 2.0 *
                                                      dot(intersect_data.hitpoint.normal, current_ray.dir));
                        // レイの屈折方向がオブジェクトの内側方向か外側方向かを確認する
                        const bool is_into = dot(intersect_data.hitpoint.normal, orienting_normal) >
                                             0.0;  // レイがオブジェクトから出るのか、入るのか

                        // Snellの法則
                        const double nc = 1.0;  // 真空の屈折率
                        // edupt::kIor オブジェクトの屈折率
                        const double nnt = is_into ? nc / edupt::kIor : edupt::kIor / nc;
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
                        const double a = edupt::kIor - nc, b = edupt::kIor + nc;
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

        std::cout << "Done. (" << photon_map->Size() << " photons are stored)" << std::endl;
        std::cout << "Creating KD-tree..." << std::endl;
        photon_map->CreateKDtree();
        std::cout << "Done." << std::endl;
    }

    int render(const std::string& filename, const int width, const int height, const int samples,
               const int supersamples, int photon_num = 5000, double gather_photon_radius = 32.0,
               int gahter_max_photon_num = 64)
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
        // スクリーンを張るベクトル
        const Vec screen_x = normalize(cross(camera_dir, camera_up)) * screen_width;
        const Vec screen_y = normalize(cross(screen_x, camera_dir)) * screen_height;
        const Vec screen_center = camera_position + camera_dir * screen_dist;

        Color* image = new Color[width * height];

        std::cout << width << "x" << height << " " << samples * (supersamples * supersamples) << " spp" << std::endl;

        PhotonMap photon_map;
        create_photonmap(photon_num, &photon_map, spheres[LightID]);

        ValueSampler<double> sampler01(0, 1);
        // OpenMP
        // // #pragma omp parallel for schedule(dynamic, 1) num_threads(4)

        for (int y = 0; y < height; y++)
        {
            std::cerr << "Rendering (y = " << y << ") " << (100.0 * y / (height - 1)) << "%" << std::endl;

            for (int x = 0; x < width; x++)
            {
                const int image_index = (height - y - 1) * width + x;
                // supersamples x supersamples のスーパーサンプリング
                for (int sy = 0; sy < supersamples; sy++)
                {
                    for (int sx = 0; sx < supersamples; sx++)
                    {
                        Color accumulated_radiance = Color();
                        // 一つのサブピクセルあたりsamples回サンプリングする
                        for (int s = 0; s < samples; s++)
                        {
                            const double rate = (1.0 / supersamples);
                            const double r1 = sx * rate + rate / 2.0;
                            const double r2 = sy * rate + rate / 2.0;
                            // スクリーン上の位置
                            const Vec screen_position = screen_center + screen_x * ((r1 + x) / width - 0.5) +
                                                        screen_y * ((r2 + y) / height - 0.5);
                            // レイを飛ばす方向
                            const Vec dir = normalize(screen_position - camera_position);

                            accumulated_radiance =
                                accumulated_radiance + photonmap_radiance(Ray(camera_position, dir), sampler01, 0,
                                                                          &photon_map, gather_photon_radius,
                                                                          gahter_max_photon_num) /
                                                           samples / (supersamples * supersamples);
                        }
                        image[image_index] = image[image_index] + accumulated_radiance;
                    }
                }
            }
        }

        //   #pragma omp parallel for schedule(dynamic, 1)

        // for (int y = 0; y < height; y++)
        // {
        //     std::cerr << "Rendering " << (100.0 * y / (height - 1)) << "%" << std::endl;
        //     srand(y * y * y);
        //     for (int x = 0; x < width; x++)
        //     {
        //         int image_index = y * width + x;
        //         image[image_index] = Color();

        //         // 2x2のサブピクセルサンプリング
        //         for (int sy = 0; sy < 2; sy++)
        //         {
        //             for (int sx = 0; sx < 2; sx++)
        //             {
        //                 // テントフィルターによってサンプリング
        //                 //
        //                 //
        //                 ピクセル範囲で一様にサンプリングするのではなく、ピクセル中央付近にサンプルがたくさん集まるように偏りを生じさせる
        //                 // const double
        //                 const double r1 = 2.0 * sampler01.sample(),
        //                              dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
        //                 const double r2 = 2.0 * sampler01.sample(),
        //                              dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);
        //                 Vec dir = screen_x * (((sx + 0.5 + dx) / 2.0 + x) / width - 0.5) +
        //                           screen_y * (((sy + 0.5 + dy) / 2.0 + y) / height - 0.5) + camera_dir;
        //                 image[image_index] = image[image_index] +
        //                                      photonmap_radiance(Ray(camera_position, dir), sampler01, 0, &photon_map,
        //                                                         gather_photon_radius, gahter_max_photon_num);
        //             }
        //         }
        //     }
        // }


        // 出力
        //        save_ppm_file(std::string("image.ppm"), image, width, height);
        save_hdr_file(filename, image, width, height);
        return 0;
    }
}  // namespace photonmap