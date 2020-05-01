#pragma once

// Copyright (c) 2014 hole
// This software is released under the MIT License (http://kagamin.net/hole/license.txt).
// A part of this software is based on smallpt (http://www.kevinbeason.com/smallpt/) and
// released under the MIT License (http://kagamin.net/hole/smallpt-license.txt).
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "intersection.h"
#include "sampler/value_sampler.h"
#include "scene.h"

namespace photonmap_sample
{
    const double PI = 3.14159265358979323846;
    const double INF = 1e20;
    const double EPS = 1e-6;
    const double MaxDepth = 5;

    // *** その他の関数 ***
    inline double rand01() { return (double)rand() / RAND_MAX; }

    // *** データ構造 ***
    struct Vec
    {
        double x, y, z;
        Vec(const double x_ = 0, const double y_ = 0, const double z_ = 0) : x(x_), y(y_), z(z_) {}
        inline Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
        inline Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
        inline Vec operator*(const double b) const { return Vec(x * b, y * b, z * b); }
        inline Vec operator/(const double b) const { return Vec(x / b, y / b, z / b); }
        inline const double length_squared() const { return x * x + y * y + z * z; }
        inline const double Length() const { return sqrt(length_squared()); }
    };
    inline Vec operator*(double f, const Vec &v) { return v * f; }
    inline Vec Normalize(const Vec &v) { return v / v.Length(); }
    // 要素ごとの積をとる
    inline const Vec Multiply(const Vec &v1, const Vec &v2) { return Vec(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z); }
    inline const double Dot(const Vec &v1, const Vec &v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }
    inline const Vec Cross(const Vec &v1, const Vec &v2)
    {
        return Vec((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
    }
    typedef Vec Color;
    const Color BackgroundColor(0.0, 0.0, 0.0);

    struct Ray
    {
        Vec org, dir;
        Ray(const Vec org_, const Vec &dir_) : org(org_), dir(dir_) {}
    };

    enum ReflectionType
    {
        DIFFUSE,     // 完全拡散面。いわゆるLambertian面。
        SPECULAR,    // 理想的な鏡面。
        REFRACTION,  // 理想的なガラス的物質。
    };

    struct Sphere
    {
        double radius;
        Vec position;
        Color emission, color;
        ReflectionType ref_type;

        Sphere(const double radius_, const Vec &position_, const Color &emission_, const Color &color_,
               const ReflectionType ref_type_)
            : radius(radius_), position(position_), emission(emission_), color(color_), ref_type(ref_type_)
        {
        }
        // 入力のrayに対する交差点までの距離を返す。交差しなかったら0を返す。
        const double intersect(const Ray &ray)
        {
            Vec o_p = position - ray.org;
            const double b = Dot(o_p, ray.dir), det = b * b - Dot(o_p, o_p) + radius * radius;
            if (det >= 0.0)
            {
                const double sqrt_det = sqrt(det);
                const double t1 = b - sqrt_det, t2 = b + sqrt_det;
                if (t1 > EPS)
                    return t1;
                else if (t2 > EPS)
                    return t2;
            }
            return 0.0;
        }
    };

    // 以下フォトンマップ用データ構造
    struct Photon
    {
        Vec position;
        Color power;
        Vec incident;

        Photon(const Vec &position_, const Color &power_, const Vec &incident_)
            : position(position_), power(power_), incident(incident_)
        {
        }
    };

    // KD-tree
    template <typename T>
    class KDTree
    {
    public:
        // k-NN searchのクエリ
        struct Query
        {
            double max_distance2;   // 探索の最大半径
            size_t max_search_num;  // 最大探索点数
            Vec search_position;    // 探索中心
            Vec normal;             // 探索中心における法線
            Query(const Vec &search_position_, const Vec &normal_, const double max_distance2_,
                  const size_t max_search_num_)
                : max_distance2(max_distance2_),
                  normal(normal_),
                  max_search_num(max_search_num_),
                  search_position(search_position_)
            {
            }
        };
        // 結果のQueueに乗せるためのデータ構造。
        struct ElementForQueue
        {
            const T *point;
            double distance2;
            ElementForQueue(const T *point_, const double distance2_) : point(point_), distance2(distance2_) {}
            bool operator<(const ElementForQueue &b) const { return distance2 < b.distance2; }
        };
        // KNNの結果を格納するキュー
        typedef std::priority_queue<ElementForQueue, std::vector<ElementForQueue> > ResultQueue;

    private:
        std::vector<T> points;
        struct KDTreeNode
        {
            T *point;
            KDTreeNode *left;
            KDTreeNode *right;
            int axis;
        };
        KDTreeNode *root;
        void delete_kdtree(KDTreeNode *node)
        {
            if (node == NULL) return;
            delete_kdtree(node->left);
            delete_kdtree(node->right);
            delete node;
        }

        // フツーのk-NN search。
        void locate_points(typename KDTree<T>::ResultQueue *pqueue, KDTreeNode *node, typename KDTree<T>::Query &query)
        {
            if (node == NULL) return;
            const int axis = node->axis;

            double delta;
            switch (axis)
            {
                case 0:
                    delta = query.search_position.x - node->point->position.x;
                    break;
                case 1:
                    delta = query.search_position.y - node->point->position.y;
                    break;
                case 2:
                    delta = query.search_position.z - node->point->position.z;
                    break;
            }

            // 対象点<->探索中心の距離が設定半径以下　かつ　対象点<->探索中心の法線方向の距離が一定以下　という条件ならその対象点格納
            const Vec dir = node->point->position - query.search_position;
            const double distance2 = dir.length_squared();
            const double dt = Dot(query.normal, dir / sqrt(distance2));
            if (distance2 < query.max_distance2 && fabs(dt) <= query.max_distance2 * 0.01)
            {
                pqueue->push(ElementForQueue(node->point, distance2));
                if (pqueue->size() > query.max_search_num)
                {
                    pqueue->pop();
                    query.max_distance2 = pqueue->top().distance2;
                }
            }
            if (delta > 0.0)
            {  // みぎ
                locate_points(pqueue, node->right, query);
                if (delta * delta < query.max_distance2)
                {
                    locate_points(pqueue, node->left, query);
                }
            }
            else
            {  // ひだり
                locate_points(pqueue, node->left, query);
                if (delta * delta < query.max_distance2)
                {
                    locate_points(pqueue, node->right, query);
                }
            }
        }

        static bool kdtree_less_operator_x(const T &left, const T &right) { return left.position.x < right.position.x; }
        static bool kdtree_less_operator_y(const T &left, const T &right) { return left.position.y < right.position.y; }
        static bool kdtree_less_operator_z(const T &left, const T &right) { return left.position.z < right.position.z; }

        KDTreeNode *create_kdtree_sub(typename std::vector<T>::iterator begin, typename std::vector<T>::iterator end,
                                      int depth)
        {
            if (end - begin <= 0)
            {
                return NULL;
            }
            const int axis = depth % 3;
            // 中央値
            switch (axis)
            {
                case 0:
                    std::sort(begin, end, kdtree_less_operator_x);
                    break;
                case 1:
                    std::sort(begin, end, kdtree_less_operator_y);
                    break;
                case 2:
                    std::sort(begin, end, kdtree_less_operator_z);
                    break;
            }
            const int median = (end - begin) / 2;
            KDTreeNode *node = new KDTreeNode;
            node->axis = axis;
            node->point = &(*(begin + median));
            // 子供
            node->left = create_kdtree_sub(begin, begin + median, depth + 1);
            node->right = create_kdtree_sub(begin + median + 1, end, depth + 1);
            return node;
        }

    public:
        KDTree() { root = NULL; }
        virtual ~KDTree() { delete_kdtree(root); }
        size_t Size() { return points.size(); }
        void SearchKNN(typename KDTree::ResultQueue *pqueue, typename KDTree<T>::Query &query)
        {
            locate_points(pqueue, root, query);
        }
        void AddPoint(const T &point) { points.push_back(point); }
        void CreateKDtree() { root = create_kdtree_sub(points.begin(), points.end(), 0); }
    };
    typedef KDTree<Photon> PhotonMap;

    // *** レンダリングするシーンデータ ****
    // from smallpt
    // Sphere spheres[] = {
    //     Sphere(5.0, Vec(50.0, 75.0, 81.6), Color(12, 12, 12), Color(), DIFFUSE),             //照明
    //     Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Color(), Color(0.75, 0.25, 0.25), DIFFUSE),    // 左
    //     Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Color(), Color(0.25, 0.25, 0.75), DIFFUSE),  // 右
    //     Sphere(1e5, Vec(50, 40.8, 1e5), Color(), Color(0.75, 0.75, 0.75), DIFFUSE),          // 奥
    //     Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Color(), Color(), DIFFUSE),                   // 手前
    //     Sphere(1e5, Vec(50, 1e5, 81.6), Color(), Color(0.75, 0.75, 0.75), DIFFUSE),          // 床
    //     Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Color(), Color(0.75, 0.75, 0.75), DIFFUSE),  // 天井
    //     Sphere(16.5, Vec(27, 16.5, 47), Color(), Color(1, 1, 1) * .99, SPECULAR),            // 鏡
    //     Sphere(16.5, Vec(73, 16.5, 78), Color(), Color(1, 1, 1) * .99, REFRACTION),          //ガラス
    // };

    Sphere spheres[] = {
        Sphere(5, Vec(50.0, 75.0, 81.6), Color(36, 36, 36), Color(), DIFFUSE),               //照明
        Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Color(), Color(0.75, 0.25, 0.25), DIFFUSE),    // 左
        Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Color(), Color(0.25, 0.25, 0.75), DIFFUSE),  // 右
        Sphere(1e5, Vec(50, 40.8, 1e5), Color(), Color(0.75, 0.75, 0.75), DIFFUSE),          // 奥
        Sphere(1e5, Vec(50, 40.8, -1e5 + 250), Color(), Color(), DIFFUSE),                   // 手前
        Sphere(1e5, Vec(50, 1e5, 81.6), Color(), Color(0.75, 0.75, 0.75), DIFFUSE),          // 床
        Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Color(), Color(0.75, 0.75, 0.75), DIFFUSE),  // 天井
        Sphere(20, Vec(65, 20, 20), Color(), Color(0.25, 0.75, 0.25), DIFFUSE),              // 緑球
        Sphere(16.5, Vec(27, 16.5, 47), Color(), Color(0.99, 0.99, 0.99), SPECULAR),         // 鏡
        Sphere(16.5, Vec(77, 16.5, 78), Color(), Color(0.99, 0.99, 0.99), REFRACTION),       //ガラス
    };

    const int LightID = 0;

    // *** レンダリング用関数 ***
    // シーンとの交差判定関数
    inline bool intersect_scene(const Ray &ray, double *t, int *id)
    {
        const double n = sizeof(spheres) / sizeof(Sphere);
        *t = INF;
        *id = -1;
        for (int i = 0; i < int(n); i++)
        {
            double d = spheres[i].intersect(ray);
            if (d > 0.0 && d < *t)
            {
                *t = d;
                *id = i;
            }
        }
        return *t < INF;
    }

    // フォトン追跡法によりフォトンマップ構築
    void create_photon_map(const int shoot_photon_num, PhotonMap *photon_map)
    {
        std::cout << "Shooting photons... (" << shoot_photon_num << " photons)" << std::endl;
        for (int i = 0; i < shoot_photon_num; i++)
        {
            // 光源からフォトンを発射する
            // 光源上の一点をサンプリングする
            const double r1 = 2 * M_PI * rand01();
            const double r2 = 1.0 - 2.0 * rand01();
            const Vec light_pos =
                spheres[LightID].position + ((spheres[LightID].radius + EPS) *
                                             Vec(sqrt(1.0 - r2 * r2) * cos(r1), sqrt(1.0 - r2 * r2) * sin(r1), r2));

            const Vec normal = Normalize(light_pos - spheres[LightID].position);
            // 光源上の点から半球サンプリングする
            Vec w, u, v;
            w = normal;
            if (fabs(w.x) > 0.1)
                u = Normalize(Cross(Vec(0.0, 1.0, 0.0), w));
            else
                u = Normalize(Cross(Vec(1.0, 0.0, 0.0), w));
            v = Cross(w, u);
            // コサイン項に比例させる。フォトンが運ぶのが放射輝度ではなく放射束であるため。
            const double u1 = 2 * PI * rand01();
            const double u2 = rand01(), u2s = sqrt(u2);
            Vec light_dir = Normalize((u * cos(u1) * u2s + v * sin(u1) * u2s + w * sqrt(1.0 - u2)));

            Ray now_ray(light_pos, light_dir);
            // emissionの値は放射輝度だが、フォトンが運ぶのは放射束なので変換する必要がある。
            // L（放射輝度）= dΦ/(cosθdωdA)なので、光源の放射束はΦ =
            // ∫∫L・cosθdωdAになる。今回は球光源で完全拡散光源であることから
            // 球上の任意の場所、任意の方向に等しい放射輝度Leを持つ。（これがemissionの値）よって、
            // Φ =
            // Le・∫∫cosθdωdAで、Le・∫dA∫cosθdωとなり、∫dAは球の面積なので4πr^2、∫cosθdωは立体角の積分なのでπとなる。
            // よって、Φ =
            // Le・4πr^2・πとなる。この値を光源から発射するフォトン数で割ってやれば一つのフォトンが運ぶ放射束が求まる。
            Color now_flux =
                spheres[LightID].emission * 4.0 * PI * pow(spheres[LightID].radius, 2.0) * PI / shoot_photon_num;

            // フォトンがシーンを飛ぶ
            bool trace_end = false;
            for (; !trace_end;)
            {
                // 放射束が0.0なフォトンを追跡してもしょうがないので打ち切る
                if (std::max(now_flux.x, std::max(now_flux.y, now_flux.z)) <= 0.0) break;

                double t;  // レイからシーンの交差位置までの距離
                int id;    // 交差したシーン内オブジェクトのID
                if (!intersect_scene(now_ray, &t, &id)) break;
                const Sphere &obj = spheres[id];
                const Vec hitpoint = now_ray.org + t * now_ray.dir;     // 交差位置
                const Vec normal = Normalize(hitpoint - obj.position);  // 交差位置の法線
                const Vec orienting_normal = Dot(normal, now_ray.dir) < 0.0
                                                 ? normal
                                                 : (-1.0 * normal);  // 交差位置の法線（物体からのレイの入出を考慮）

                switch (obj.ref_type)
                {
                    case DIFFUSE:
                    {
                        // 拡散面なのでフォトンをフォトンマップに格納する
                        photon_map->AddPoint(Photon(hitpoint, now_flux, now_ray.dir));

                        // 反射するかどうかをロシアンルーレットで決める
                        // 例によって確率は任意。今回はフォトンマップ本に従ってRGBの反射率の平均を使う
                        const double probability = (obj.color.x + obj.color.y + obj.color.z) / 3;
                        if (probability > rand01())
                        {  // 反射
                            // orienting_normalの方向を基準とした正規直交基底(w, u,
                            // v)を作る。この基底に対する半球内で次のレイを飛ばす。
                            Vec w, u, v;
                            w = orienting_normal;
                            if (fabs(w.x) > 0.1)
                                u = Normalize(Cross(Vec(0.0, 1.0, 0.0), w));
                            else
                                u = Normalize(Cross(Vec(1.0, 0.0, 0.0), w));
                            v = Cross(w, u);
                            // コサイン項を使った重点的サンプリング
                            const double r1 = 2 * PI * rand01();
                            const double r2 = rand01(), r2s = sqrt(r2);
                            Vec dir = Normalize((u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.0 - r2)));

                            now_ray = Ray(hitpoint, dir);
                            now_flux = Multiply(now_flux, obj.color) / probability;
                            continue;
                        }
                        else
                        {  // 吸収（すなわちここで追跡終了）
                            trace_end = true;
                            continue;
                        }
                    }
                    break;
                    case SPECULAR:
                    {
                        // 完全鏡面なのでフォトン格納しない
                        // 完全鏡面なのでレイの反射方向は決定的。
                        now_ray = Ray(hitpoint, now_ray.dir - normal * 2.0 * Dot(normal, now_ray.dir));
                        now_flux = Multiply(now_flux, obj.color);
                        continue;
                    }
                    break;
                    case REFRACTION:
                    {
                        // やはりフォトン格納しない
                        Ray reflection_ray = Ray(hitpoint, now_ray.dir - normal * 2.0 * Dot(normal, now_ray.dir));
                        bool into = Dot(normal, orienting_normal) > 0.0;  // レイがオブジェクトから出るのか、入るのか

                        // Snellの法則
                        const double nc = 1.0;  // 真空の屈折率
                        const double nt = 1.5;  // オブジェクトの屈折率
                        const double nnt = into ? nc / nt : nt / nc;
                        const double ddn = Dot(now_ray.dir, orienting_normal);
                        const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

                        if (cos2t < 0.0)
                        {  // 全反射した
                            now_ray = reflection_ray;
                            now_flux = Multiply(now_flux, obj.color);
                            continue;
                        }
                        // 屈折していく方向
                        Vec tdir =
                            Normalize(now_ray.dir * nnt - normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)));

                        // SchlickによるFresnelの反射係数の近似
                        const double a = nt - nc, b = nt + nc;
                        const double R0 = (a * a) / (b * b);
                        const double c = 1.0 - (into ? -ddn : Dot(tdir, normal));
                        const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
                        const double Tr = 1.0 - Re;  // 屈折光の運ぶ光の量
                        const double probability = Re;

                        // 屈折と反射のどちらか一方を追跡する。
                        // ロシアンルーレットで決定する。
                        if (rand01() < probability)
                        {  // 反射
                            now_ray = reflection_ray;
                            // Fresnel係数Reを乗算し、ロシアンルーレット確率prob.で割る。
                            // 今、prob.=Reなので Re / prob. = 1.0 となる。
                            // よって、now_flux = Multiply(now_flux, obj.color) * Re / probability; が以下の式になる。
                            // 屈折の場合も同様。
                            now_flux = Multiply(now_flux, obj.color);
                            continue;
                        }
                        else
                        {  // 屈折
                            now_ray = Ray(hitpoint, tdir);
                            now_flux = Multiply(now_flux, obj.color);
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

    // ray方向からの放射輝度を求める
    Color radiance(const Ray &ray, const int depth, PhotonMap *photon_map, const double gather_radius,
                   const int gahter_max_photon_num)
    {
        double t;  // レイからシーンの交差位置までの距離
        int id;    // 交差したシーン内オブジェクトのID
        if (!intersect_scene(ray, &t, &id)) return BackgroundColor;

        const Sphere &obj = spheres[id];
        const Vec hitpoint = ray.org + t * ray.dir;             // 交差位置
        const Vec normal = Normalize(hitpoint - obj.position);  // 交差位置の法線
        const Vec orienting_normal =
            Dot(normal, ray.dir) < 0.0 ? normal : (-1.0 * normal);  // 交差位置の法線（物体からのレイの入出を考慮）

        // 色の反射率最大のものを得る。ロシアンルーレットで使う。
        // ロシアンルーレットの閾値は任意だが色の反射率等を使うとより良い。
        double russian_roulette_probability = std::max(obj.color.x, std::max(obj.color.y, obj.color.z));
        // 一定以上レイを追跡したらロシアンルーレットを実行し追跡を打ち切るかどうかを判断する
        if (depth > MaxDepth)
        {
            if (rand01() >= russian_roulette_probability) return obj.emission;
        }
        else
            russian_roulette_probability = 1.0;  // ロシアンルーレット実行しなかった

        switch (obj.ref_type)
        {
            case DIFFUSE:
            {
                // フォトンマップをつかって放射輝度推定する
                PhotonMap::ResultQueue pqueue;
                // k近傍探索。gather_radius半径内のフォトンを最大gather_max_photon_num個集めてくる
                PhotonMap::Query query(hitpoint, orienting_normal, gather_radius, gahter_max_photon_num);
                photon_map->SearchKNN(&pqueue, query);
                Color accumulated_flux;
                double max_distance2 = -1;

                // キューからフォトンを取り出しvectorに格納する
                std::vector<PhotonMap::ElementForQueue> photons;
                photons.reserve(pqueue.size());
                for (; !pqueue.empty();)
                {
                    PhotonMap::ElementForQueue p = pqueue.top();
                    pqueue.pop();
                    photons.push_back(p);
                    max_distance2 = std::max(max_distance2, p.distance2);
                }

                // 円錐フィルタを使用して放射輝度推定する
                const double max_distance = sqrt(max_distance2);
                const double k = 1.1;
                for (int i = 0; i < photons.size(); i++)
                {
                    const double w = 1.0 - (sqrt(photons[i].distance2) / (k * max_distance));  // 円錐フィルタの重み
                    const Color v = Multiply(obj.color, photons[i].point->power) /
                                    PI;  // Diffuse面のBRDF = 1.0 / πであったのでこれをかける
                    accumulated_flux = accumulated_flux + w * v;
                }
                accumulated_flux = accumulated_flux / (1.0 - 2.0 / (3.0 * k));  // 円錐フィルタの係数
                if (max_distance2 > 0.0)
                {
                    return obj.emission + accumulated_flux / (PI * max_distance2) / russian_roulette_probability;
                }
            }
            break;

                // SPECULARとREFRACTIONの場合はパストレーシングとほとんど変わらない。
                // 単純に反射方向や屈折方向の放射輝度(Radiance)をradiance()で求めるだけ。

            case SPECULAR:
            {
                // 完全鏡面にヒットした場合、反射方向から放射輝度をもらってくる
                return obj.emission + radiance(Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir)), depth + 1,
                                               photon_map, gather_radius, gahter_max_photon_num) /
                                          russian_roulette_probability;
            }
            break;
            case REFRACTION:
            {
                Ray reflection_ray = Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir));
                bool into = Dot(normal, orienting_normal) > 0.0;  // レイがオブジェクトから出るのか、入るのか

                // Snellの法則
                const double nc = 1.0;  // 真空の屈折率
                const double nt = 1.5;  // オブジェクトの屈折率
                const double nnt = into ? nc / nt : nt / nc;
                const double ddn = Dot(ray.dir, orienting_normal);
                const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

                if (cos2t < 0.0)
                {  // 全反射した
                    // 反射方向から放射輝度をもらってくる
                    return obj.emission +
                           Multiply(obj.color, radiance(Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir)),
                                                        depth + 1, photon_map, gather_radius, gahter_max_photon_num)) /
                               russian_roulette_probability;
                }
                // 屈折していく方向
                Vec tdir = Normalize(ray.dir * nnt - normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)));

                // SchlickによるFresnelの反射係数の近似
                const double a = nt - nc, b = nt + nc;
                const double R0 = (a * a) / (b * b);
                const double c = 1.0 - (into ? -ddn : Dot(tdir, normal));
                const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
                const double Tr = 1.0 - Re;  // 屈折光の運ぶ光の量
                const double probability = 0.25 + 0.5 * Re;

                // 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
                // ロシアンルーレットで決定する。
                if (depth > 2)
                {
                    if (rand01() < probability)
                    {  // 反射
                        return obj.emission + Multiply(obj.color, radiance(reflection_ray, depth + 1, photon_map,
                                                                           gather_radius, gahter_max_photon_num) *
                                                                      Re) /
                                                  probability / russian_roulette_probability;
                    }
                    else
                    {  // 屈折
                        return obj.emission + Multiply(obj.color, radiance(Ray(hitpoint, tdir), depth + 1, photon_map,
                                                                           gather_radius, gahter_max_photon_num) *
                                                                      Tr) /
                                                  (1.0 - probability) / russian_roulette_probability;
                    }
                }
                else
                {  // 屈折と反射の両方を追跡
                    return obj.emission + Multiply(obj.color, radiance(reflection_ray, depth + 1, photon_map,
                                                                       gather_radius, gahter_max_photon_num) *
                                                                      Re +
                                                                  radiance(Ray(hitpoint, tdir), depth + 1, photon_map,
                                                                           gather_radius, gahter_max_photon_num) *
                                                                      Tr) /
                                              russian_roulette_probability;
                }
            }
            break;
        }

        return Color();
    }

    // コサイン項によるImportance Sampling
    decltype(auto) cosine_sampling(const Vec &w, const Vec &u, const Vec &v, ValueSampler<double> &sampler01)
    {
        const double u1 = 2 * M_PI * sampler01.sample();
        const double u2 = sampler01.sample();
        const double u3 = std::sqrt(u2);
        return Normalize((u * std::cos(u1) * u3 + v * std::sin(u1) * u3 + w * std::sqrt(1.0 - u2)));
    }
    // 球の１点をサンプリング
    decltype(auto) sphere_sampling(ValueSampler<double> &sampler01)
    {
        const double r1 = 2 * M_PI * sampler01.sample();
        const double r2 = 1.0 - 2.0 * sampler01.sample();
        const double r3 = 1.0 - r2 * r2;
        return Vec(std::sqrt(r3) * std::cos(r1), std::sqrt(r3) * std::sin(r1), r2);
    }

    void create_photonmap(const int photon_count, PhotonMap *photon_map, const Sphere &light_sphere)
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
            const Vec source_dir = Normalize(source_pos - light_sphere.position);

            // 光源上の点から半球サンプリングする
            Vec w, u, v;
            w = source_dir;

            if (fabs(w.x) > 0.1)
            {
                u = Normalize(Cross(Vec(0.0, 1.0, 0.0), w));
            }
            else
            {
                u = Normalize(Cross(Vec(1.0, 0.0, 0.0), w));
            }
            v = Cross(w, u);

            Vec light_dir = cosine_sampling(w, u, v, sampler01);

            Ray current_ray(source_pos, light_dir);

            Color current_flux = light_sphere.emission * 4.0 * std::pow(light_sphere.radius * M_PI, 2.0) / photon_count;

            bool trace_end = false;
            while (!trace_end)
            {
                if (std::max(current_flux.x, std::max(current_flux.y, current_flux.z)) <= 0.0) break;

                double t;  // レイからシーンの交差位置までの距離
                int id;    // 交差したシーン内オブジェクトのID
                if (!intersect_scene(current_ray, &t, &id)) break;
                const Sphere &obj = spheres[id];

                const Vec hitpoint = current_ray.org + t * current_ray.dir;  // 交差位置
                const Vec normal = Normalize(hitpoint - obj.position);       // 交差位置の法線
                const Vec orienting_normal = Dot(normal, current_ray.dir) < 0.0
                                                 ? normal
                                                 : (-1.0 * normal);  // 交差位置の法線（物体からのレイの入出を考慮）

                switch (obj.ref_type)
                {
                    case DIFFUSE:
                    {
                        // 拡散面なのでフォトンをフォトンマップに格納する
                        photon_map->AddPoint(Photon(hitpoint, current_flux, current_ray.dir));

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
                                diffuse_u = Normalize(Cross(Vec(0.0, 1.0, 0.0), diffuse_w));
                            }
                            else
                            {
                                diffuse_u = Normalize(Cross(Vec(1.0, 0.0, 0.0), diffuse_w));
                            }
                            diffuse_v = Cross(diffuse_w, diffuse_u);
                            Vec dir = cosine_sampling(diffuse_w, diffuse_u, diffuse_v, sampler01);

                            current_ray = Ray(hitpoint, dir);
                            current_flux = Multiply(current_flux, obj.color) / probability;
                        }
                        else
                        {  // 吸収
                            trace_end = true;
                        }
                    }
                    break;
                    case SPECULAR:
                    {
                        // 完全鏡面
                        current_ray = Ray(hitpoint, current_ray.dir - hitpoint * 2.0 * Dot(hitpoint, current_ray.dir));
                        current_flux = Multiply(current_flux, obj.color);
                    }
                    break;
                    case REFRACTION:
                    {
                        // 屈折
                        Ray reflection_ray =
                            Ray(hitpoint, current_ray.dir - hitpoint * 2.0 * Dot(hitpoint, current_ray.dir));
                        // レイの屈折方向がオブジェクトの内側方向か外側方向かを確認する
                        const bool is_into =
                            Dot(normal, orienting_normal) > 0.0;  // レイがオブジェクトから出るのか、入るのか

                        // Snellの法則
                        const double nc = 1.0;  // 真空の屈折率
                        // edupt::kIor オブジェクトの屈折率
                        const double nnt = is_into ? nc / edupt::kIor : edupt::kIor / nc;
                        const double ddn = Dot(current_ray.dir, orienting_normal);
                        const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

                        if (cos2t < 0.0)
                        {  // 全反射
                            current_ray = reflection_ray;
                            current_flux = Multiply(current_flux, obj.color);
                            continue;
                        }
                        // 屈折していく方向
                        Vec tdir = Normalize(current_ray.dir * nnt -
                                             normal * (is_into ? 1.0 : -1.0) * (ddn * nnt + std::sqrt(cos2t)));

                        // SchlickによるFresnelの反射係数の近似
                        const double a = edupt::kIor - nc, b = edupt::kIor + nc;
                        const double R0 = (a * a) / (b * b);
                        const double c = 1.0 - (is_into ? -ddn : Dot(tdir, normal));
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
                            current_flux = Multiply(current_flux, obj.color);
                            continue;
                        }
                        else
                        {  // 屈折
                            current_ray = Ray(hitpoint, tdir);
                            current_flux = Multiply(current_flux, obj.color);
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
    // *** .hdrフォーマットで出力するための関数 ***
    struct HDRPixel
    {
        unsigned char r, g, b, e;
        HDRPixel(const unsigned char r_ = 0, const unsigned char g_ = 0, const unsigned char b_ = 0,
                 const unsigned char e_ = 0)
            : r(r_), g(g_), b(b_), e(e_){};
        unsigned char get(int idx)
        {
            switch (idx)
            {
                case 0:
                    return r;
                case 1:
                    return g;
                case 2:
                    return b;
                case 3:
                    return e;
            }
            return 0;
        }
    };

    // doubleのRGB要素を.hdrフォーマット用に変換
    HDRPixel get_hdr_pixel(const Color &color)
    {
        double d = std::max(color.x, std::max(color.y, color.z));
        if (d <= 1e-32) return HDRPixel();
        int e;
        double m = frexp(d, &e);  // d = m * 2^e
        d = m * 256.0 / d;
        return HDRPixel(color.x * d, color.y * d, color.z * d, e + 128);
    }

    // 書き出し用関数
    void save_hdr_file(const std::string &filename, const Color *image, const int width, const int height)
    {
        FILE *fp = fopen(filename.c_str(), "wb");
        if (fp == NULL)
        {
            std::cerr << "Error: " << filename << std::endl;
            return;
        }
        // .hdrフォーマットに従ってデータを書きだす
        // ヘッダ
        unsigned char ret = 0x0a;
        fprintf(fp, "#?RADIANCE%c", (unsigned char)ret);
        fprintf(fp, "# Made with 100%% pure HDR Shop%c", ret);
        fprintf(fp, "FORMAT=32-bit_rle_rgbe%c", ret);
        fprintf(fp, "EXPOSURE=1.0000000000000%c%c", ret, ret);

        // 輝度値書き出し
        fprintf(fp, "-Y %d +X %d%c", height, width, ret);
        for (int i = height - 1; i >= 0; i--)
        {
            std::vector<HDRPixel> line;
            for (int j = 0; j < width; j++)
            {
                HDRPixel p = get_hdr_pixel(image[j + i * width]);
                line.push_back(p);
            }
            fprintf(fp, "%c%c", 0x02, 0x02);
            fprintf(fp, "%c%c", (width >> 8) & 0xFF, width & 0xFF);
            for (int i = 0; i < 4; i++)
            {
                for (int cursor = 0; cursor < width;)
                {
                    const int cursor_move = std::min(127, width - cursor);
                    fprintf(fp, "%c", cursor_move);
                    for (int j = cursor; j < cursor + cursor_move; j++) fprintf(fp, "%c", line[j].get(i));
                    cursor += cursor_move;
                }
            }
        }

        fclose(fp);
    }

    int render()
    {
        int width = 640;
        int height = 480;
        int photon_num = 500000;
        double gather_photon_radius = 32.0;
        int gahter_max_photon_num = 64;

        // カメラ位置
        Ray camera(Vec(50.0, 52.0, 295.6), Normalize(Vec(0.0, -0.042612, -1.0)));
        // シーン内でのスクリーンのx,y方向のベクトル
        Vec cx = Vec(width * 0.5135 / height);
        Vec cy = Normalize(Cross(cx, camera.dir)) * 0.5135;
        Color *image = new Color[width * height];

        // フォトンマップ構築
        PhotonMap photon_map;
        create_photon_map(photon_num, &photon_map);
        // create_photonmap(photon_num, &photon_map, spheres[LightID]);

        // // #pragma omp parallel for schedule(dynamic, 1)
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
        //                 ピクセル範囲で一様にサンプリングするのではなく、ピクセル中央付近にサンプルがたくさん集まるように偏りを生じさせる
        //                 const double r1 = 2.0 * rand01(), dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
        //                 const double r2 = 2.0 * rand01(), dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);
        //                 Vec dir = cx * (((sx + 0.5 + dx) / 2.0 + x) / width - 0.5) +
        //                           cy * (((sy + 0.5 + dy) / 2.0 + y) / height - 0.5) + camera.dir;
        //                 image[image_index] =
        //                     image[image_index] + radiance(Ray(camera.org + dir * 130.0, Normalize(dir)), 0,
        //                     &photon_map,
        //                                                   gather_photon_radius, gahter_max_photon_num);
        //             }
        //         }
        //     }
        // }
        int supersamples = 2;
        int samples = 2;
        const Vec camera_up = Vec(0.0, 1.0, 0.0);
        int screen_width = width;
        int screen_height = height;

        const Vec screen_x = Normalize(Cross(camera.dir, camera_up)) * screen_width;
        const Vec screen_y = Normalize(Cross(screen_x, camera.dir)) * screen_height;
        const Vec screen_center = camera.org + camera.dir * 800;

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
                            const Vec dir = Normalize(screen_position - camera.org);

                            accumulated_radiance = accumulated_radiance +
                                                   radiance(Ray(camera.org + dir * 130.0, Normalize(dir)), 0,
                                                            &photon_map, gather_photon_radius, gahter_max_photon_num)

                                ;
                        }
                        image[image_index] = image[image_index] + accumulated_radiance;
                    }
                }
            }
        }

        // .hdrフォーマットで出力
        save_hdr_file(std::string("image_sample.hdr"), image, width, height);
    }
}  // namespace photonmap_sample