// expanded smallppm (code is exactly the same as smallppm.cpp but with more comments)

#include <cmath>    // smallppm, Progressive Photon Mapping by T. Hachisuka
#include <cstdlib>  // originally smallpt, a path tracer by Kevin Beason, 2008
#include <fstream>
#include <iostream>  // Usage: ./smallppm 100000 && xv image.ppm

namespace photonmap
{
    namespace smallppm
    {
        constexpr double PI = 3.14159265358979;  // ^^^^^^:number of photons emitted
        constexpr double ALPHA = 0.7;            // the alpha parameter of PPM

        // Halton sequence with reverse permutation
        int primes[61] = {2,   3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,
                          59,  61,  67,  71,  73,  79,  83,  89,  97,  101, 103, 107, 109, 113, 127, 131,
                          137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
                          227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283};
        inline int rev(const int i, const int p)
        {
            if (i == 0)
                return i;
            else
                return p - i;
        }
        double hal(const int b, int j)
        {
            const int p = primes[b];
            double h = 0.0, f = 1.0 / (double)p, fct = f;
            while (j > 0)
            {
                h += rev(j % p, p) * fct;
                j /= p;
                fct *= f;
            }
            return h;
        }

        struct Vec
        {
            double x, y, z;  // vector: position, also color (r,g,b)
            Vec(double x_ = 0, double y_ = 0, double z_ = 0)
            {
                x = x_;
                y = y_;
                z = z_;
            }
            inline Vec operator+(const Vec& b) const { return Vec(x + b.x, y + b.y, z + b.z); }
            inline Vec operator-(const Vec& b) const { return Vec(x - b.x, y - b.y, z - b.z); }
            inline Vec operator+(double b) const { return Vec(x + b, y + b, z + b); }
            inline Vec operator-(double b) const { return Vec(x - b, y - b, z - b); }
            inline Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
            inline Vec mul(const Vec& b) const { return Vec(x * b.x, y * b.y, z * b.z); }
            inline Vec norm() { return (*this) * (1.0 / sqrt(x * x + y * y + z * z)); }
            inline double dot(const Vec& b) const { return x * b.x + y * b.y + z * b.z; }
            Vec operator%(Vec& b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
        };

#define MAX(x, y) ((x > y) ? x : y)

        struct AABB
        {
            Vec min, max;  // axis aligned bounding box
            inline void fit(const Vec& p)
            {
                if (p.x < min.x) min.x = p.x;  // min
                if (p.y < min.y) min.y = p.y;  // min
                if (p.z < min.z) min.z = p.z;  // min
                max.x = MAX(p.x, max.x);
                max.y = MAX(p.y, max.y);
                max.z = MAX(p.z, max.z);
            }
            inline void reset()
            {
                min = Vec(1e20, 1e20, 1e20);
                max = Vec(-1e20, -1e20, -1e20);
            }
        };

        struct HPoint
        {
            Vec f, pos, nrm, flux;
            double r2;
            unsigned int n;  // n = N / ALPHA in the paper
            int pix;
        };

        struct List
        {
            HPoint* id;
            List* next;
        };
        List* ListAdd(HPoint* i, List* h)
        {
            List* p = new List;
            p->id = i;
            p->next = h;
            return p;
        }

        unsigned int num_hash, pixel_index, num_photon;
        double hash_s;
        List** hash_grid;
        List* hitpoints = NULL;
        AABB hpbbox;

        // spatial hash function
        inline unsigned int hash(const int ix, const int iy, const int iz)
        {
            return (unsigned int)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % num_hash;
        }

        void build_hash_grid(const int w, const int h)
        {
            // find the bounding box of all the measurement points
            hpbbox.reset();
            List* lst = hitpoints;
            while (lst != NULL)
            {
                HPoint* hp = lst->id;
                lst = lst->next;
                hpbbox.fit(hp->pos);
            }

            // heuristic for initial radius
            Vec ssize = hpbbox.max - hpbbox.min;
            double irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;
            std::cout << irad << std::endl;

            // determine hash table size
            // we now find the bounding box of all the measurement points inflated by the initial radius
            hpbbox.reset();
            lst = hitpoints;
            int vphoton = 0;
            while (lst != NULL)
            {
                HPoint* hp = lst->id;
                lst = lst->next;
                hp->r2 = irad * irad;
                hp->n = 0;
                hp->flux = Vec();
                vphoton++;
                hpbbox.fit(hp->pos - irad);
                hpbbox.fit(hp->pos + irad);
            }

            // make each grid cell two times larger than the initial radius
            hash_s = 1.0 / (irad * 2.0);
            num_hash = vphoton;

            // build the hash table
            hash_grid = new List*[num_hash];
            for (unsigned int i = 0; i < num_hash; i++) hash_grid[i] = NULL;
            lst = hitpoints;
            while (lst != NULL)
            {
                HPoint* hp = lst->id;
                lst = lst->next;
                Vec BMin = ((hp->pos - irad) - hpbbox.min) * hash_s;
                Vec BMax = ((hp->pos + irad) - hpbbox.min) * hash_s;
                for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++)
                {
                    for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++)
                    {
                        for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++)
                        {
                            int hv = hash(ix, iy, iz);
                            hash_grid[hv] = ListAdd(hp, hash_grid[hv]);
                        }
                    }
                }
            }
        }

        struct Ray
        {
            Vec o, d;
            Ray(){};
            Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
        };

        enum Refl_t
        {
            DIFF,
            SPEC,
            REFR
        };  // material types, used in radiance()
        struct Sphere
        {
            double rad;
            Vec p, c;
            Refl_t refl;
            Sphere(double r_, Vec p_, Vec c_, Refl_t re_) : rad(r_), p(p_), c(c_), refl(re_) {}
            inline double intersect(const Ray& r) const
            {
                // ray-sphere intersection returns distance
                Vec op = p - r.o;
                double t, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
                if (det < 0)
                {
                    return 1e20;
                }
                else
                {
                    det = sqrt(det);
                }
                return (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : 1e20);
            }
        };

        Sphere sph[] = {  // Scene: radius, position, color, material
            Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(.75, .25, .25), DIFF),    // Left
            Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(.25, .25, .75), DIFF),  // Right
            Sphere(1e5, Vec(50, 40.8, 1e5), Vec(.75, .75, .75), DIFF),          // Back
            Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), DIFF),                // Front
            Sphere(1e5, Vec(50, 1e5, 81.6), Vec(.75, .75, .75), DIFF),          // Bottomm
            Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(.75, .75, .75), DIFF),  // Top
            Sphere(16.5, Vec(27, 16.5, 47), Vec(1, 1, 1) * .999, SPEC),         // Mirror
            Sphere(16.5, Vec(73, 16.5, 88), Vec(1, 1, 1) * .999, REFR),         // Glass
            Sphere(8.5, Vec(50, 8.5, 60), Vec(1, 1, 1) * .999, DIFF)};          // Middle
        // tone mapping and gamma correction
        int toInt(double x) { return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5); }

        // find the closet interection
        inline bool intersect(const Ray& r, double& t, int& id)
        {
            int n = sizeof(sph) / sizeof(Sphere);
            double d, inf = 1e20;
            t = inf;
            for (int i = 0; i < n; i++)
            {
                d = sph[i].intersect(r);
                if (d < t)
                {
                    t = d;
                    id = i;
                }
            }
            return t < inf;
        }

        // generate a photon ray from the point light source with QMC
        void genp(Ray* pr, Vec* f, int i)
        {
            *f = Vec(2500, 2500, 2500) * (PI * 4.0);  // flux
            double p = 2. * PI * hal(0, i), t = 2. * acos(sqrt(1. - hal(1, i)));
            double st = sin(t);
            pr->d = Vec(cos(p) * st, cos(t), sin(p) * st);
            pr->o = Vec(50, 60, 85);
        }

        void trace(const Ray& r, int dpt, bool m, const Vec& fl, const Vec& adj, int i)
        {
            double t;
            int id;

            dpt++;
            if (!intersect(r, t, id) || (dpt >= 20)) return;

            int d3 = dpt * 3;
            const Sphere& obj = sph[id];
            Vec x = r.o + r.d * t, n = (x - obj.p).norm(), f = obj.c;
            Vec nl = n.dot(r.d) < 0 ? n : n * -1;
            double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;

            if (obj.refl == DIFF)
            {
                // Lambertian

                // use QMC to sample the next direction
                double r1 = 2. * PI * hal(d3 - 1, i), r2 = hal(d3 + 0, i);
                double r2s = sqrt(r2);
                Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
                Vec v = w % u, d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

                if (m)
                {
                    // eye ray
                    // store the measurment point
                    HPoint* hp = new HPoint;
                    hp->f = f.mul(adj);
                    hp->pos = x;
                    hp->nrm = n;
                    hp->pix = pixel_index;
                    hitpoints = ListAdd(hp, hitpoints);
                }
                else
                {
                    // photon ray
                    // find neighboring measurement points and accumulate flux via progressive density estimation
                    Vec hh = (x - hpbbox.min) * hash_s;
                    int ix = abs(int(hh.x)), iy = abs(int(hh.y)), iz = abs(int(hh.z));
                    // strictly speaking, we should use #pragma omp critical here.
                    // it usually works without an artifact due to the fact that photons are
                    // rarely accumulated to the same measurement points at the same time (especially with QMC).
                    // it is also significantly faster.
                    {
                        List* hp = hash_grid[hash(ix, iy, iz)];
                        while (hp != NULL)
                        {
                            HPoint* hitpoint = hp->id;
                            hp = hp->next;
                            Vec v = hitpoint->pos - x;
                            // check normals to be closer than 90 degree (avoids some edge brightning)
                            if ((hitpoint->nrm.dot(n) > 1e-3) && (v.dot(v) <= hitpoint->r2))
                            {
                                // unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
                                double g = (hitpoint->n * ALPHA + ALPHA) / (hitpoint->n * ALPHA + 1.0);
                                hitpoint->r2 = hitpoint->r2 * g;
                                hitpoint->n++;
                                hitpoint->flux = (hitpoint->flux + hitpoint->f.mul(fl) * (1. / PI)) * g;
                            }
                        }
                    }
                    if (hal(d3 + 1, i) < p) trace(Ray(x, d), dpt, m, f.mul(fl) * (1. / p), adj, i);
                }
            }
            else if (obj.refl == SPEC)
            {
                // mirror
                trace(Ray(x, r.d - n * 2.0 * n.dot(r.d)), dpt, m, f.mul(fl), f.mul(adj), i);
            }
            else
            {
                // glass
                Ray lr(x, r.d - n * 2.0 * n.dot(r.d));
                bool into = (n.dot(nl) > 0.0);
                double nc = 1.0, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;

                // total internal reflection
                if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) return trace(lr, dpt, m, fl, adj, i);

                Vec td = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
                double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : td.dot(n));
                double Re = R0 + (1 - R0) * c * c * c * c * c, P = Re;
                Ray rr(x, td);
                Vec fa = f.mul(adj);
                if (m)
                {
                    // eye ray (trace both rays)
                    trace(lr, dpt, m, fl, fa * Re, i);
                    trace(rr, dpt, m, fl, fa * (1.0 - Re), i);
                }
                else
                {
                    // photon ray (pick one via Russian roulette)
                    (hal(d3 - 1, i) < P) ? trace(lr, dpt, m, fl, fa, i) : trace(rr, dpt, m, fl, fa, i);
                }
            }
        }

        int render()
        {
            // samps * 1000 photon paths will be traced
            int w = 1024, h = 768, samps = 1000;

            // trace eye rays and store measurement points
            Ray cam(Vec(50, 48, 295.6), Vec(0, -0.042612, -1).norm());
            Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, *c = new Vec[w * h], vw;
            for (int y = 0; y < h; y++)
            {
                fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0 * y / (h - 1));
                for (int x = 0; x < w; x++)
                {
                    pixel_index = x + y * w;
                    Vec d = cx * ((x + 0.5) / w - 0.5) + cy * (-(y + 0.5) / h + 0.5) + cam.d;
                    trace(Ray(cam.o + d * 140, d.norm()), 0, true, Vec(), Vec(1, 1, 1), 0);
                }
            }
            fprintf(stderr, "\n");

            // build the hash table over the measurement points
            build_hash_grid(w, h);

            // trace photon rays with multi-threading
            num_photon = samps;
            vw = Vec(1, 1, 1);
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < num_photon; i++)
            {
                double p = 100. * (i + 1) / num_photon;
                fprintf(stderr, "\rPhotonPass %5.2f%%", p);
                int m = 1000 * i;
                Ray r;
                Vec f;
                for (int j = 0; j < 1000; j++)
                {
                    genp(&r, &f, m + j);
                    trace(r, 0, 0 > 1, f, vw, m + j);
                }
            }

            // density estimation
            List* lst = hitpoints;
            while (lst != NULL)
            {
                HPoint* hp = lst->id;
                lst = lst->next;
                int i = hp->pix;
                c[i] = c[i] + hp->flux * (1.0 / (PI * hp->r2 * num_photon * 1000.0));
            }

            // save the image after tone mapping and gamma correction
            FILE* f = fopen("image.ppm", "w");
            fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
            for (int i = 0; i < w * h; i++)
            {
                fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
            }
        }
    }  // namespace smallppm
}  // namespace photonmap