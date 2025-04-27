#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#include <memory>

#include "texture.h"
#include "color.h"
#include "rtw_stb_image.h"

class texture {
    public:
        virtual ~texture() = default;

        virtual color value(double u, double v, const point3& p) const = 0;
};

class solid_color : public texture {
    public:
        solid_color(const color& albedo) : albedo(albedo) {}

        solid_color(double red, double green, double blue) : solid_color(color(red,green,blue)) {}

        color value(double u, double v, const point3& p) const override {
            return albedo;
        }

    private:
        color albedo;
};

class checker_texture : public texture {
    public:
        checker_texture(double scale, std::shared_ptr<texture> even, std::shared_ptr<texture> odd) 
            : inv_scale(1.0 / scale), even(even), odd(odd) {}
        checker_texture(double scale, const color& c1, const color& c2)
            : checker_texture(scale, std::make_shared<solid_color>(c1), std::make_shared<solid_color>(c2)) {}
        color value(double u, double v, const point3& p) const override {
            // return color(0.5, 0.5, 0.5);
            auto xInt = int(std::floor(inv_scale * p.x()));
            auto yInt = int(std::floor(inv_scale * p.y()));
            auto zInt = int(std::floor(inv_scale * p.z()));
            
            bool isEven = (xInt + yInt + zInt) % 2 == 0;

            return isEven ? even->value(u, v, p) : odd->value(u, v, p);
        }
    
    private:
        double inv_scale;
        std::shared_ptr<texture> even;
        std::shared_ptr<texture> odd;
};

class image_texture : public texture {
    public:
        image_texture(const char* filename) : image(filename) {}

        color value(double u, double v, const point3& p) const override {
            if (image.height() <= 0) return color(0, 1, 1); // Cyan texture if there is no texture data available

            /* Make sure to clamp the u and v between 0 and 1 */
            u = interval(0,1).clamp(u);
            v = 1.0 - interval(0,1).clamp(v);

            /* Get the pixel value at some point in the image texture */
            auto i = int(u * image.width());
            auto j = int(v * image.height());
            auto pixel = image.pixel_data(i,j);

            auto color_scale = 1.0 / 255.0;
            return color(color_scale*pixel[0], color_scale*pixel[1], color_scale*pixel[2]);
        }
    private:
        rtw_image image;
};

#endif