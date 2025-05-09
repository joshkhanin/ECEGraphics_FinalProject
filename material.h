#ifndef MATERIAL_H
#define MATERIAL_H
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
#include "texture.h"

class material {
  public:
    virtual ~material() = default;

    virtual bool scatter(
        const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
    ) const {
        return false;
    }

    virtual color emitted(double u, double v, const point3& p) const {
      return color(0,0,0);
    }
};


class lambertian : public material {
  public:
    // lambertian(const color& albedo) : albedo(albedo) {}
    lambertian(const color& albedo) : tex(make_shared<solid_color>(albedo)) {}
    lambertian(shared_ptr<texture> tex) : tex(tex) {}

    bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered)
    const override {
        auto scatter_direction = rec.normal + random_unit_vector();

        // Catch degenerate scatter direction
        if (scatter_direction.near_zero())
            scatter_direction = rec.normal;

        scattered = ray(rec.p, scatter_direction);
        // attenuation = albedo;
        attenuation = tex->value(rec.u, rec.v, rec.p);
        return true;
    }

  private:
    shared_ptr<texture> tex;
};


class metal : public material {
  public:
    metal(const color& albedo, double fuzz) : albedo(albedo), fuzz(fuzz < 1 ? fuzz : 1) {}

    bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered)
    const override {
        vec3 reflected = reflect(r_in.direction(), rec.normal);
        reflected = unit_vector(reflected) + (fuzz * random_unit_vector());
        scattered = ray(rec.p, reflected);
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
    }

  private:
    color albedo;
    double fuzz;
};


class dielectric : public material {
  public:
    dielectric(double refraction_index) : refraction_index(refraction_index) {}

    bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered)
    const override {
        attenuation = color(1.0, 1.0, 1.0);
        double ri = rec.front_face ? (1.0/refraction_index) : refraction_index;

        vec3 unit_direction = unit_vector(r_in.direction());
        double cos_theta = std::fmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);

        bool cannot_refract = ri * sin_theta > 1.0;
        vec3 direction;

        if (cannot_refract || reflectance(cos_theta, ri) > random_double())
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, ri);

        scattered = ray(rec.p, direction);
        return true;
    }

  private:
    // Refractive index in vacuum or air, or the ratio of the material's refractive index over
    // the refractive index of the enclosing media
    double refraction_index;

    static double reflectance(double cosine, double refraction_index) {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1 - refraction_index) / (1 + refraction_index);
        r0 = r0*r0;
        return r0 + (1-r0)*std::pow((1 - cosine),5);
    }
};

class diffuse_light : public material {
  public:
    diffuse_light(shared_ptr<texture> tex, double intensity) : tex(tex), intensity(intensity) {}
    color emitted(double u, double v, const point3& p) const override {
          return tex->value(u, v, p) * intensity;
    }
  private:
    shared_ptr<texture> tex;
    double intensity;
};

class earth_mat : public material {
  public:
    earth_mat(shared_ptr<texture> albedo, shared_ptr<texture> specular) : albedo(albedo), specular(specular) {}

    bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered)
    const override {
        double reflect_amount = specular->value(rec.u, rec.v, rec.p).length_squared() / 9.0;
        double sample_threshold = reflect_amount / 2;

        if (random_double() < sample_threshold) {
          /* Create a normal vector perturbed according to the nature of ocean waves */
          auto gaussian_vec = random_gaussian_vector(0.5);
          auto squished_vec = gaussian_vec - dot(gaussian_vec, rec.normal) * gaussian_vec;
          auto perturbed_normal = unit_vector(squished_vec + rec.normal);
          auto scatter_direction = unit_vector(reflect(r_in.direction(), perturbed_normal));
          scattered = ray(rec.p, scatter_direction);
          // attenuation = albedo;
          attenuation = color(1.0, 1.0, 1.0) * 0.3;
          return (dot(scattered.direction(), rec.normal) > 0);
        } else {
          /* Reflect according to a Lambertian */
          auto scatter_direction = rec.normal + random_unit_vector();

          // Catch degenerate scatter direction
          if (scatter_direction.near_zero())
              scatter_direction = rec.normal;

          scattered = ray(rec.p, scatter_direction);
          // attenuation = albedo;
          attenuation = albedo->value(rec.u, rec.v, rec.p) * (1.0 / ((1-sample_threshold)+0.05));
          return true;
        }
        
        // TODO: Adjust scatter direction based on albedo texture
        // auto scatter_direction = perturbed_normal + random_unit_vector();
        // auto scatter_direction = unit_vector(reflect(r_in.direction(), perturbed_normal));


        // // Catch degenerate scatter direction
        // if (scatter_direction.near_zero())
        //     scatter_direction = rec.normal;

        // scattered = ray(rec.p, scatter_direction);
        // // attenuation = albedo;
        // attenuation = albedo->value(rec.u, rec.v, rec.p);
        // return (dot(scattered.direction(), rec.normal) > 0);
        // return true;
    }
  private:
    shared_ptr<texture> albedo;
    shared_ptr<texture> specular;
};

// TODO: Add an Earth material, which will store multiple textures:
//    - Albedo
//    - Reflectivity
//    - Night lighting


#endif
