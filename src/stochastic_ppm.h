#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

#include "cameras.h"
#include "hdr.h"
#include "intersection.h"
#include "kdtree.h"
#include "material.h"
#include "ppm.h"
#include "radiance.h"
#include "random.h"
#include "ray.h"
#include "sampler/value_sampler.h"
#include "scene.h"
#include "vec.h"

namespace photonmap
{
    namespace sppm
    {
        using namespace edupt;
        using namespace photonmap::utility;
        template <int Width, int Height>
        class StochasticPpm
        {
        public:
            int LightID = 0;
            const double INF = 1e20;
            const double EPS = 1e-6;
            int min_depth = 5;
            int depth_threshold = 24;
            double alpha = 0.7;            // the alpha parameter of PPM
            double initial_radius_2 = 25;  // the alpha parameter of PPM
            int* imageidx_pointidx;

            // カメラ位置
            const Vec camera_position = Vec(50.0, 52.0, 220.0);
            const Vec camera_dir = normalize(Vec(0.0, -0.04, -1.0));
            const Vec camera_up = Vec(0.0, 1.0, 0.0);

            // ワールド座標系でのスクリーンの大きさ
            double screen_width;
            double screen_height = 30.0;
            // スクリーンまでの距離
            double screen_dist = 40.0;

            const int width = Width;
            const int height = Height;
            static constexpr int PixelCount = Width * Height;

            int samples;
            int supersamples;
            int photon_num;

            int valid_photon = 0;

            std::array<Color, PixelCount> ld_container;
            std::array<double, PixelCount> phi_container;
            std::array<Color, PixelCount> tau_container;

            std::random_device seed_gen_;
            std::mt19937 engine_ = std::mt19937(seed_gen_());

            /*
            DoFCamera camera;
            DoFCamera InitCamera()
            {
                return DoFCamera(width, height, screen_height, screen_dist, camera_position, camera_dir, camera_up,
                                 supersamples, 2, 60, engine_);
            }
            /*/
            PinholeCamera camera;
            PinholeCamera InitCamera()
            {
                return PinholeCamera(width, height, screen_height, screen_dist, camera_position, camera_dir, camera_up,
                                     supersamples);
            }

            // */
            StochasticPpm() : screen_width(30.0 * width / height), screen_height(30.0), camera(InitCamera()) {}
            StochasticPpm(const int samples_, const int supersamples_)
                : samples(samples_),
                  supersamples(supersamples_),
                  screen_width(30.0 * width / height),
                  screen_height(30.0),
                  camera(InitCamera())
            {
            }
            struct ProgressiveIntersection
            {
                Color emission;
                // SHARED
                Color accumulated_flux;
                int index;
                int photon_count;
                double photon_radius_2;  // photon_radius^2

                Intersection intersection;
                Vec position;
                double weight;
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

            struct SPPMPixel
            {
                // SHARED
                Color accumulated_flux;
                int index;
                double photon_count;
                int photon_add_count;
                double photon_radius_2;  // photon_radius^2

                bool is_visible = false;

                Intersection intersection;
                Vec position;
                Color weight;
                SPPMPixel(Intersection intersection_, const Color& weight_, double photon_radius_2_,
                          double photon_count_, int index_)
                    : intersection(intersection_),
                      position(intersection_.hitpoint.position),
                      weight(weight_),
                      photon_radius_2(photon_radius_2_),
                      index(index_)
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

            using PointMap = KDTree<SPPMPixel>;

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

            void reset_ld_container()
            {
                for (size_t i = 0; i < PixelCount; i++)
                {
                    ld_container[i] = Color();
                }
            }

            // ray方向からの放射輝度を求める
            void create_point_loop(const Ray& input_ray, ValueSampler<double>* rnd, PointMap* point_map, int index)
            {
                const int LightID = 0;
                const double INF = 1e20;
                const double EPS = 1e-6;
                const double MaxDepth = 5;
                int depth = -1;
                bool direct = true;
                Color weight(1, 1, 1);
                Color rad(0, 0, 0);

                Ray ray(input_ray.org, input_ray.dir);

                while (true)
                {
                    depth++;
                    Intersection intersection;
                    // シーンと交差判定
                    if (!intersect_scene(ray, &intersection))
                    {
                        break;
                    }

                    const Sphere& now_object = spheres[intersection.object_id];
                    const Hitpoint& hitpoint = intersection.hitpoint;
                    const Vec orienting_normal =
                        dot(hitpoint.normal, ray.dir) < 0.0
                            ? hitpoint.normal
                            : (-1.0 * hitpoint.normal);  // 交差位置の法線（物体からのレイの入出を考慮）

                    if (direct)
                    {
                        ld_container[index] = ld_container[index] + multiply(weight, now_object.emission);
                    }

                    direct = now_object.reflection_type != REFLECTION_TYPE_DIFFUSE;

                    switch (now_object.reflection_type)
                    {
                        // 完全拡散面
                        case REFLECTION_TYPE_DIFFUSE:
                        {
                            auto est = next_event_est(intersection, spheres[LightID], orienting_normal, *rnd);

                            ld_container[index] = ld_container[index] + multiply(weight, est);
                            auto p_it = imageidx_pointidx[index];

                            // Not Found
                            if (p_it < 0)
                            {
                                std::cout << point_map->Size() << std::endl;
                                imageidx_pointidx[index] = 1;
                                point_map->AddData(SPPMPixel(intersection, weight, initial_radius_2, 0, index));
                                break;
                            }
                            else  // Found Point
                            {
                                auto& point_list = point_map->GetData();
                                for (auto& point : point_list)
                                {
                                    if (point.index == index)
                                    {
                                        point.intersection = intersection;
                                        point.position = intersection.hitpoint.position;
                                        point.weight = weight;
                                        // ld_container[index] = ld_container[index] +
                                        //                         now_object.emission * weight /
                                        //                         russian_roulette_probability;
                                        break;
                                    }
                                }
                            }

                            // // orienting_normalの方向を基準とした正規直交基底(w, u,
                            // // v)を作る。この基底に対する半球内で次のレイを飛ばす。
                            // Vec w, u, v;
                            // w = orienting_normal;
                            // if (fabs(w.x) >
                            //     kEPS)  //
                            //     ベクトルwと直交するベクトルを作る。w.xが0に近い場合とそうでない場合とで使うベクトルを変える。
                            //     u = normalize(cross(Vec(0.0, 1.0, 0.0), w));
                            // else
                            //     u = normalize(cross(Vec(1.0, 0.0, 0.0), w));
                            // v = cross(w, u);
                            // // コサイン項を使った重点的サンプリング
                            // const double r1 = 2 * kPI * rnd->sample();
                            // const double r2 = rnd->sample(), r2s = sqrt(r2);
                            // Vec dir = normalize((u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.0 - r2)));

                            // ray = Ray(hitpoint.position, dir);
                            weight = multiply(weight, now_object.color);

                            return;
                        }
                        break;

                        // 完全鏡面
                        case REFLECTION_TYPE_SPECULAR:
                        {
                            ray =
                                Ray(hitpoint.position, ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir));
                            weight = multiply(weight, now_object.color);
                        }
                        break;

                        // 屈折率kIorのガラス
                        case REFLECTION_TYPE_REFRACTION:
                        {
                            const Ray reflection_ray =
                                Ray(hitpoint.position, ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir));
                            const bool into = dot(hitpoint.normal, orienting_normal) >
                                              0.0;  // レイがオブジェクトから出るのか、入るのか

                            // Snellの法則
                            const double nc = 1.0;   // 真空の屈折率
                            const double nt = kIor;  // オブジェクトの屈折率
                            const double nnt = into ? nc / nt : nt / nc;
                            const double ddn = dot(ray.dir, orienting_normal);
                            const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

                            if (cos2t < 0.0)
                            {  // 全反射
                                ray = reflection_ray;
                                weight = multiply(weight, now_object.color);
                                break;
                            }

                            // 屈折の方向
                            const Ray refraction_ray =
                                Ray(hitpoint.position, normalize(ray.dir * nnt - hitpoint.normal * (into ? 1.0 : -1.0) *
                                                                                     (ddn * nnt + sqrt(cos2t))));

                            // SchlickによるFresnelの反射係数の近似を使う
                            const double a = nt - nc, b = nt + nc;
                            const double R0 = (a * a) / (b * b);

                            const double c = 1.0 - (into ? -ddn : dot(refraction_ray.dir, -1.0 * orienting_normal));
                            const double Re =
                                R0 +
                                (1.0 - R0) *
                                    pow(c,
                                        5.0);  // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
                            const double nnt2 = pow(
                                into ? nc / nt : nt / nc,
                                2.0);  // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
                            const double Tr = (1.0 - Re) * nnt2;  // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

                            // 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
                            // ロシアンルーレットで決定する。
                            const double probability = 0.25 + 0.5 * Re;

                            if (rnd->sample() < probability)
                            {  // 反射
                                ray = reflection_ray;
                                weight = multiply(weight, now_object.color * Re / (probability));
                            }
                            else
                            {
                                ray = refraction_ray;
                                weight = multiply(weight, now_object.color * Tr / ((1.0 - probability)));
                            }
                        }
                        break;
                    }
                    // 色の反射率最大のものを得る。ロシアンルーレットで使う。
                    // ロシアンルーレットの閾値は任意だが色の反射率等を使うとより良い。
                    double russian_roulette_probability =
                        std::max(now_object.color.x, std::max(now_object.color.y, now_object.color.z));

                    // 反射回数が一定以上になったらロシアンルーレットの確率を急上昇させる。（スタックオーバーフロー対策）
                    if (depth >= MaxDepth) russian_roulette_probability *= std::pow(0.5, depth - MaxDepth);

                    // ロシアンルーレットを実行し追跡を打ち切るかどうかを判断する。
                    // ただしDepth回の追跡は保障する。
                    if (depth >= kDepth)
                    {
                        if (rnd->sample() >= russian_roulette_probability) break;
                    }
                    else
                        russian_roulette_probability = 1.0;  // ロシアンルーレット実行しなかった

                    weight = weight / russian_roulette_probability;
                }
            }

            bool create_photon(const Sphere& light_sphere, ValueSampler<double>& sampler01, PointMap* point_map,
                               double gather_radius, double gahter_max_photon_num)
            {
                // 光源からフォトンを発射する
                // 光源は球。球の一点をサンプリングする

                const Vec source_pos =
                    light_sphere.position + ((light_sphere.radius + EPS) * sphere_sampling(sampler01));
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

                Color weight = light_sphere.emission * 4.0 * std::pow(light_sphere.radius * M_PI, 2.0);

                bool trace_end = false;
                bool added_photon = false;
                int depth = -1;

                while (!trace_end)
                {
                    if (std::max(weight.x, std::max(weight.y, weight.z)) <= 0.0) break;

                    Intersection intersect_data;
                    if (!intersect_scene(current_ray, &intersect_data)) break;
                    const Sphere& obj = spheres[intersect_data.object_id];
                    // 物体に対する入射方向が表か裏かを確認する
                    const Vec orienting_normal = dot(intersect_data.hitpoint.normal, current_ray.dir) < 0.0
                                                     ? intersect_data.hitpoint.normal  // 物体外から入射
                                                     : (-1.0 * intersect_data.hitpoint.normal);  // 物体中から入射
                    depth++;
                    switch (obj.reflection_type)
                    {
                        case edupt::REFLECTION_TYPE_DIFFUSE:
                        {
                            if (depth > 0)  // add contribution (depth == 0 => direct illumination)
                            {
                                typename PointMap::ResultQueue result_queue;

                                typename PointMap::Query query(intersect_data.hitpoint.position, orienting_normal,
                                                               gather_radius, gahter_max_photon_num);
                                point_map->SearchNearest(&result_queue, query);
                                // キューからフォトンを取り出しvectorに格納する
                                std::vector<typename PointMap::ElementForQueue> points;
                                points.reserve(result_queue.size());

                                while (!result_queue.empty())
                                {
                                    typename PointMap::ElementForQueue p = result_queue.top();
                                    result_queue.pop();
                                    points.push_back(p);
                                }

                                for (auto& point : points)
                                {
                                    auto diff = point.point->position - intersect_data.hitpoint.position;
                                    if (dot(point.point->intersection.hitpoint.normal,
                                            intersect_data.hitpoint.normal) <= 1e-3 ||
                                        dot(diff, diff) > point.point->photon_radius_2)
                                        continue;

                                    // Update _pixel_ $\Phi$ and $M$ for nearby
                                    // photon
                                    auto f = spheres[point.point->intersection.object_id].color /
                                             M_PI;  // point.point->color / M_PI /(M_PI / )

                                    Color accumulated_flux = multiply(weight, f);
                                    point.point->accumulated_flux = point.point->accumulated_flux + accumulated_flux;
                                    point.point->photon_add_count++;
                                    added_photon = true;
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

                                const double r1 = 2 * M_PI * sampler01.sample();
                                const double r2 = sampler01.sample(), r2s = std::sqrt(r2);
                                Vec dir = normalize(
                                    (u * std::cos(r1) * r2s + v * std::sin(r1) * r2s + w * std::sqrt(1.0 - r2)));

                                auto& f_pdf_cos = obj.color;  // cos * point.point->color / M_PI /(cos / M_PI);
                                auto fl_new = multiply(
                                    weight, f_pdf_cos);  // * std::abs(dot(dir, intersect_data.hitpoint.normal));

                                weight = fl_new / probability;

                                current_ray = Ray(intersect_data.hitpoint.position, dir);
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
                            current_ray =
                                Ray(intersect_data.hitpoint.position,
                                    current_ray.dir - intersect_data.hitpoint.normal * 2.0 *
                                                          dot(intersect_data.hitpoint.normal, current_ray.dir));
                            weight = multiply(weight, obj.color);
                            continue;
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
                                weight = multiply(weight, obj.color);
                                continue;
                            }
                            // 屈折していく方向
                            Vec tdir = normalize(current_ray.dir * nnt - intersect_data.hitpoint.normal *
                                                                             (is_into ? 1.0 : -1.0) *
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
                                // よって、weight = Multiply(weight, obj.color) * Re / probability;
                                // が以下の式になる。 屈折の場合も同様。
                                weight = multiply(weight, obj.color);
                                continue;
                            }
                            else
                            {  // 屈折
                                current_ray = Ray(intersect_data.hitpoint.position, tdir);
                                weight = multiply(weight, obj.color);
                                continue;
                            }
                        }

                        break;
                    }
                }

                return added_photon;
            }
            void update_photons_with_render(int photon_count, PointMap* point_map, double gather_radius,
                                            double gahter_max_photon_num, int width, int height,
                                            ValueSampler<double>& sampler01)
            {
                double max_radius_2;
                for (const auto& p : point_map->GetData())
                {
                    max_radius_2 = std::max(max_radius_2, p.photon_radius_2);
                }
                std::cout << "Rmax2 = " << max_radius_2 << std::endl;

                for (size_t i = 0; i < photon_count; i++)
                {
                    create_photon(spheres[LightID], sampler01, point_map, max_radius_2, gahter_max_photon_num);
                }
            };

            int render(const std::string& filename, const int width, const int height, const int samples,
                       const int supersamples, int photon_num, double gather_photon_radius, int gahter_max_photon_num,
                       int eye_pass_count = 100)
            {
                // // カメラ位置
                // const Vec camera_position = Vec(50.0, 52.0, 220.0);
                // const Vec camera_dir = normalize(Vec(0.0, -0.04, -1.0));
                // const Vec camera_up = Vec(0.0, 1.0, 0.0);

                // // ワールド座標系でのスクリーンの大きさ
                // const double screen_width = 30.0 * width / height;
                // const double screen_height = 30.0;
                // // スクリーンまでの距離
                // const double screen_dist = 40.0;
                // // スクリーンを張るベクトル
                // const Vec screen_x = normalize(cross(camera_dir, camera_up)) * screen_width;
                // const Vec screen_y = normalize(cross(screen_x, camera_dir)) * screen_height;
                // const Vec screen_center = camera_position + camera_dir * screen_dist;

                Color* image = new Color[PixelCount];
                imageidx_pointidx = new int[PixelCount];
                for (size_t i = 0; i < PixelCount; i++)
                {
                    imageidx_pointidx[i] = -1;
                }
                std::vector<ProgressiveIntersection> hitpoint_list;

                std::cout << width << "x" << height << " " << samples << " spp" << std::endl;

                PointMap point_map;

                // for (size_t i = 0; i < PixelCount; i++)
                // {
                //     tau_container[i] = Color();
                // }

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
                sampler01.sample();

                std::cout << "Creating Point..." << std::endl;
                //    reset_ld_container();

                for (int i = 0; i < eye_pass_count; i++)
                {
                    for (int y = 0; y < height; y++)
                    {
                        // std::cout << "Creating Point Tracing Map (y = " << y << ") " << (100.0 * y / (height - 1))
                        //           << "%" << std::endl;

                        for (int x = 0; x < width; x++)
                        {
                            const int image_index = (height - y - 1) * width + x;

                            // for (int sy = 0; sy < supersamples; sy++)
                            // {
                            //     for (int sx = 0; sx < supersamples; sx++)
                            //     {
                            //         // 一つのサブピクセルあたりsamples回サンプリングする
                            //         for (int s = 0; s < samples; s++)
                            //         {
                            //             create_point(camera.get_ray(x, y, sx, sy), sampler01, 0, 1, &point_map,
                            //                          image_index, false);
                            //         }
                            //     }
                            // }

                            // 一つのサブピクセルあたりsamples回サンプリングする

                            int sx = sampler01.sample() * supersamples;
                            int sy = sampler01.sample() * supersamples;
                            for (int s = 0; s < samples; s++)
                            {
                                // create_point(camera.get_ray(x, y, sx, sy), sampler01, 0, 1, &point_map, image_index,
                                //              false);
                                create_point_loop(camera.get_ray(x, y, sx, sy), &sampler01, &point_map, image_index);
                            }
                        }
                    }

                    std::cout << "Creating KDTree..." << std::endl;
                    point_map.CreateKDtree();
                    std::cout << "Done Creating KDTree" << std::endl;

                    update_photons_with_render(photon_num, &point_map, gather_photon_radius, gahter_max_photon_num,
                                               width, height, sampler01);
                    //            update_photons(photon_num, &point_map, gather_photon_radius, gahter_max_photon_num);
                    std::cout << "Done Emitting Photon" << std::endl;

                    std::cout << "Rendering..." << std::endl;

                    // for (auto node : point_map.GetData())
                    // {
                    //     // image[node.index] =
                    //     //     node.emission + node.weight * node.accumulated_flux *
                    //     //                         (1.0 / (M_PI * node.photon_radius_2 * ((i + 1) * photon_num)));
                    //     image[node.index] =
                    //         image[node.index] + multiply(node.weight, node.accumulated_flux) *
                    //                                 (1.0 / (M_PI * node.photon_radius_2 * ((i + 1) * photon_num)));
                    //     if (node.index == image_index_sample)
                    //     {
                    //         std::cout << image[node.index] << std::endl;
                    //     }
                    // }

                    auto Np = (i + 1) * photon_num;
                    auto& nodes = point_map.GetData();
                    for (auto& node : nodes)
                    {
                        int p_it = node.index;
                        if (node.photon_add_count > 0)
                        {
                            double Nnew = node.photon_count + alpha * node.photon_add_count;
                            double Rnew2 = node.photon_radius_2 * Nnew / (node.photon_count + node.photon_add_count);

                            tau_container[p_it] = (tau_container[p_it] + multiply(node.accumulated_flux, node.weight)) *
                                                  Rnew2 / node.photon_radius_2;
                            node.photon_count = Nnew;
                            node.photon_radius_2 = Rnew2;
                            //                            valid_photon += node.photon_add_count;
                            node.photon_add_count = 0;
                            node.accumulated_flux = Color();
                        }

                        node.weight = Color();
                    }

                    for (size_t p_it = 0; p_it < PixelCount; p_it++)
                    {
                        image[p_it] = ld_container[p_it] / ((i + 1));
                    }
                    const int image_index_sample = 49196;

                    std::cout << image_index_sample << std::endl;

                    // for (size_t p_it = 0; p_it < PixelCount; p_it++)
                    // {
                    //     auto& node = nodes[p_it];
                    //     auto Np = (i + 1) * photon_num;
                    //     image[node.index] =
                    //         image[node.index] + tau_container[node.index] * (1.0 / (M_PI * node.photon_radius_2 *
                    //         Np));
                    // }

                    for (auto& node : nodes)
                    {
                        int p_it = node.index;
                        // auto Np = valid_photon;
                        image[p_it] = image[p_it] + tau_container[p_it] * (1.0 / (M_PI * node.photon_radius_2 * Np));
                    }

                    std::cout << image[image_index_sample] << std::endl;
                    //                    image[image_index_sample] = Color(1, 0, 0);

                    save_ppm_file(filename + "_" + std::to_string(i) + ".ppm", image, width, height);
                    save_hdr_file(filename + "_" + std::to_string(i) + ".hdr", image, width, height);

                    std::cout << "Done Rendering" << std::endl;
                }

                // 出力
                //        save_ppm_file(std::string("image.ppm"), image, width, height);
                save_ppm_file(filename + ".ppm", image, width, height);
                save_hdr_file(filename + ".hdr", image, width, height);
                return 0;
            }
        };
    }  // namespace sppm
}  // namespace photonmap