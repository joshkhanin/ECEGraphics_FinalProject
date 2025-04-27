#ifndef SPHERE_H
#define SPHERE_H
//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "hittable.h"


class sphere : public hittable {
  public:
    sphere(const point3& center, double radius, shared_ptr<material> mat)
      : center(center), radius(std::fmax(0,radius)), mat(mat) {}

    bool hit(const ray& r, interval ray_t, hit_record& rec) const override {
        vec3 oc = center - r.origin();
        auto a = r.direction().length_squared();
        auto h = dot(r.direction(), oc);
        auto c = oc.length_squared() - radius*radius;

        auto discriminant = h*h - a*c;
        if (discriminant < 0)
            return false;

        auto sqrtd = std::sqrt(discriminant);

        // Find the nearest root that lies in the acceptable range.
        auto root = (h - sqrtd) / a;
        if (!ray_t.surrounds(root)) {
            root = (h + sqrtd) / a;
            if (!ray_t.surrounds(root))
                return false;
        }

        rec.t = root;
        rec.p = r.at(rec.t);
        vec3 outward_normal = (rec.p - center) / radius;
        rec.set_face_normal(r, outward_normal);
        get_sphere_uv(outward_normal, rec.u, rec.v);
        rec.mat = mat;

        return true;
    }

    // TODO: Add function static void get_sphere_uv(const point3& p, double& u, double& v);
    static void get_sphere_uv(const point3& p, double& u, double& v) {
      /**
       * p is a point on the unit sphere (in Cartesian coordinates) and get_sphere_uv maps p to spherical angles theta and phi
       * y = cos(theta), x = sin(theta)cos(phi), z = sin(theta)cos(phi)
       * y = -cos(theta), x = -sin(theta)cos(phi), z= sin(theta)sin(phi)
       */

       auto theta = std::acos(-p.y());            v = theta / pi;
       auto phi = std::atan2(-p.z(), p.x()) + pi; u = phi / (2*pi);
    }

  private:
    point3 center;
    double radius;
    shared_ptr<material> mat;
};


#endif
