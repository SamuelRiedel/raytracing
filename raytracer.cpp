// [header]
// A very basic raytracer example.
// [/header]
// [compile]
// c++ -o raytracer -O3 -Wall raytracer.cpp
// [/compile]
// [ignore]
// Copyright (C) 2012  www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// [/ignore]
#ifdef __cplusplus
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#else
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#define M_PI 3.141592653589793
#endif

#define VEC3_xyz(X, Y, Z) ((Vec3) { .x = X, .y = Y, .z = Z })
#define VEC3_x(X) VEC3_xyz(X, X, X)
#define VEC3 VEC3_x(0)
typedef struct Vec3 {
  float x, y, z;
} Vec3;

float vec_dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

float vec_length2(Vec3 v) { return vec_dot(v, v); }
float vec_length(Vec3 v) { return sqrt(vec_length(v)); }

Vec3 vec_scale(Vec3 v, float s) {
  v.x *= s;
  v.y *= s;
  v.z *= s;
  return v;
}

Vec3 vec_add(Vec3 a, Vec3 b) {
  return (Vec3) { .x = a.x + b.x, .y = a.y + b.y, .z = a.z + b.z };
}
Vec3 vec_subtract(Vec3 a, Vec3 b) {
  return (Vec3) { .x = a.x - b.x, .y = a.y - b.y, .z = a.z - b.z };
}
Vec3 vec_negate(Vec3 v) {
  return (Vec3) { .x = -v.x, .y = -v.y, .z = -v.z };
}
Vec3 vec_multiply(Vec3 a, Vec3 b) {
  return (Vec3) { .x = a.x * b.x, .y = a.y * b.y, .z = a.z * b.z };
}

Vec3 vec_normalize(Vec3 v) {
  float nor2 = vec_length2(v);
  if (nor2 > 0) {
    float invNor = 1 / sqrt(nor2);
    v = vec_scale(v, invNor);
  }
  return v;
}

typedef struct Sphere {
  Vec3 center;                      /// position of the sphere
  float radius, radius2;            /// sphere radius and radius^2
  Vec3 surfaceColor, emissionColor; /// surface color and emission (light)
  float transparency, reflection;   /// surface transparency and reflectivity
} Sphere;

// Compute a ray-sphere intersection using the geometric solution
bool intersect(const Sphere sphere, const Vec3 rayorig, const Vec3 raydir,
               float *t0, float *t1) {
  Vec3 l = vec_subtract(sphere.center, rayorig);
  float tca = vec_dot(l, raydir);
  if (tca < 0)
    return false;
  float d2 = vec_length2(l) - tca * tca;
  if (d2 > sphere.radius2)
    return false;
  float thc = sqrt(sphere.radius2 - d2);
  *t0 = tca - thc;
  *t1 = tca + thc;

  return true;
}

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
#define MAX_RAY_DEPTH 5

float mix(const float a, const float b, const float mix) {
  return b * mix + a * (1 - mix);
}

//[comment]
// This is the main trace function. It takes a ray as argument (defined by its
// origin
// and direction). We test if this ray intersects any of the geometry in the
// scene.
// If the ray intersects an object, we compute the intersection point, the
// normal
// at the intersection point, and shade this point using this information.
// Shading depends on the surface property (is it transparent, reflective,
// diffuse).
// The function returns a color for the ray. If the ray intersects an object
// that
// is the color of the object at the intersection point, otherwise it returns
// the background color.
//[/comment]
Vec3 trace(const Vec3 rayorig, const Vec3 raydir, const Sphere *spheres,
           const unsigned num_spheres, const int depth) {
  // if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
  float tnear = INFINITY;
  const Sphere *sphere = NULL;
  // find intersection of this ray with the sphere in the scene
  for (unsigned i = 0; i < num_spheres; ++i) {
    float t0 = INFINITY, t1 = INFINITY;
    if (intersect(spheres[i], rayorig, raydir, &t0, &t1)) {
      if (t0 < 0)
        t0 = t1;
      if (t0 < tnear) {
        tnear = t0;
        sphere = &spheres[i];
      }
    }
  }
  // if there's no intersection return black or background color
  if (!sphere)
    return VEC3_x(2);
  Vec3 surfaceColor =
      VEC3; // color of the ray/surfaceof the object intersected by the ray
  Vec3 phit =
      vec_add(rayorig, vec_scale(raydir, tnear)); // point of intersection
  Vec3 nhit =
      vec_subtract(phit, sphere->center); // normal at the intersection point
  nhit = vec_normalize(nhit);             // normalize normal direction
  // If the normal and the view direction are not opposite to each other
  // reverse the normal direction. That also means we are inside the sphere so
  // set
  // the inside bool to true. Finally reverse the sign of IdotN which we want
  // positive.
  float bias = 1e-4; // add some bias to the point from which we will be tracing
  bool inside = false;
  if (vec_dot(raydir, nhit) > 0) {
    nhit = vec_negate(nhit);
    inside = true;
  }
  if ((sphere->transparency > 0 || sphere->reflection > 0) &&
      depth < MAX_RAY_DEPTH) {
    float facingratio = -vec_dot(raydir, nhit);
    // change the mix value to tweak the effect
    float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
    // compute reflection direction (not need to normalize because all vectors
    // are already normalized)
    Vec3 refldir =
        vec_subtract(raydir, vec_scale(nhit, 2 * vec_dot(raydir, nhit)));
    refldir = vec_normalize(refldir);
    Vec3 reflection = trace(vec_add(phit, vec_scale(nhit, bias)), refldir,
                            spheres, num_spheres, depth + 1);
    Vec3 refraction = VEC3;
    // if the sphere is also transparent compute refraction ray (transmission)
    if (sphere->transparency) {
      float ior = 1.1,
            eta = (inside) ? ior
                           : 1 / ior; // are we inside or outside the surface?
      float cosi = -vec_dot(nhit, raydir);
      float k = 1 - eta * eta * (1 - cosi * cosi);
      Vec3 refrdir = vec_add(vec_scale(raydir, eta),
                             vec_scale(nhit, (eta * cosi - sqrt(k))));
      refrdir = vec_normalize(refrdir);
      refraction = trace(vec_subtract(phit, vec_scale(nhit, bias)), refrdir,
                         spheres, num_spheres, depth + 1);
    }
    // the result is a mix of reflection and refraction (if the sphere is
    // transparent)
    surfaceColor =
        vec_multiply(vec_add(vec_scale(reflection, fresneleffect),
                             vec_scale(refraction, (1 - fresneleffect) *
                                                       sphere->transparency)),
                     sphere->surfaceColor);
  } else {
    // it's a diffuse object, no need to raytrace any further
    for (unsigned i = 0; i < num_spheres; ++i) {
      if (spheres[i].emissionColor.x > 0) {
        // this is a light
        float transmission = 1;
        Vec3 lightDirection = vec_subtract(spheres[i].center, phit);
        lightDirection = vec_normalize(lightDirection);
        for (unsigned j = 0; j < num_spheres; ++j) {
          if (i != j) {
            float t0, t1;
            if (intersect(spheres[j], vec_add(phit, vec_scale(nhit, bias)),
                          lightDirection, &t0, &t1)) {
              transmission = 0;
              break;
            }
          }
        }
        float dot = vec_dot(nhit, lightDirection);
        if (dot > 0) {
          transmission *= dot;
          Vec3 tmp = vec_scale(sphere->surfaceColor, transmission);
          tmp = vec_multiply(tmp, spheres[i].emissionColor);
          surfaceColor = vec_add(tmp, surfaceColor);
        }
      }
    }
  }

  return vec_add(surfaceColor, sphere->emissionColor);
}

//[comment]
// Main rendering function. We compute a camera ray for each pixel of the image
// trace it and return a color. If the ray hits a sphere, we return the color of
// the
// sphere at the intersection point, else we return the background color.
//[/comment]
void render(const Sphere *spheres, unsigned num_spheres) {
  unsigned width = 640, height = 480;
  // unsigned width = 3840, height = 2160;
  // unsigned width = 1280, height = 720;
  Vec3 image[width * height];
  Vec3 *pixel = image;
  float invWidth = 1 / (float)width;
  float invHeight = 1 / (float)height;
  float fov = 30, aspectratio = width / (float)height;
  float angle = tan(M_PI * 0.5 * fov / 180.);
  // Trace rays
  for (unsigned y = 0; y < height; ++y) {
    for (unsigned x = 0; x < width; ++x, ++pixel) {
      float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
      float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
      Vec3 raydir = VEC3_xyz(xx, yy, -1);
      raydir = vec_normalize(raydir);
      *pixel = trace(VEC3_x(0), raydir, spheres, num_spheres, 0);
    }
  }
  // Save result to a PPM image (keep these flags if you compile under Windows)
  FILE *fptr;
  fptr = fopen("./untitled.ppm", "wb");
  fprintf(fptr, "P6\n%d %d\n255\n", width, height);
  for (unsigned i = 0; i < width * height; ++i) {
    unsigned char x = image[i].x > 1 ? 255 : (unsigned char)(image[i].x * 255);
    unsigned char y = image[i].y > 1 ? 255 : (unsigned char)(image[i].y * 255);
    unsigned char z = image[i].z > 1 ? 255 : (unsigned char)(image[i].z * 255);
    fprintf(fptr, "%c%c%c", x, y, z);
  }
  fclose(fptr);
}

//[comment]
// In the main function, we will create the scene which is composed of 5 spheres
// and 1 light (which is also a sphere). Then, once the scene description is
// complete
// we render that scene, by calling the render() function.
//[/comment]
int main(int argc, char **argv) {
  const unsigned num_spheres = 6;
  Sphere spheres[6];
  // position, radius, surface color, reflectivity, transparency, emission color
  spheres[0] = (Sphere) { .center = VEC3_xyz(0.0, -10004, -20),
                          .radius = 10000,
                          .radius2 = 10000 * 10000,
                          .surfaceColor = VEC3_xyz(0.20, 0.20, 0.20),
                          .emissionColor = VEC3,
                          .transparency = 0,
                          .reflection = 0 };
  spheres[1] = (Sphere) { .center = VEC3_xyz(0.0, 0, -20),
                          .radius = 4,
                          .radius2 = 4 * 4,
                          .surfaceColor = VEC3_xyz(1.00, 0.32, 0.36),
                          .emissionColor = VEC3,
                          .transparency = 0,
                          .reflection = 1 };
  spheres[2] = (Sphere) { .center = VEC3_xyz(5.0, -1, -15),
                          .radius = 2,
                          .radius2 = 2 * 2,
                          .surfaceColor = VEC3_xyz(0.90, 0.76, 0.46),
                          .emissionColor = VEC3,
                          .transparency = 0,
                          .reflection = 0 };
  spheres[3] = (Sphere) { .center = VEC3_xyz(5.0, 0, -25),
                          .radius = 3,
                          .radius2 = 3 * 3,
                          .surfaceColor = VEC3_xyz(0.65, 0.77, 0.97),
                          .emissionColor = VEC3,
                          .transparency = 0,
                          .reflection = 1 };
  spheres[4] = (Sphere) { .center = VEC3_xyz(-5.5, 0, -15),
                          .radius = 3,
                          .radius2 = 3 * 3,
                          .surfaceColor = VEC3_xyz(0.90, 0.90, 0.90),
                          .emissionColor = VEC3,
                          .transparency = 0,
                          .reflection = 1 };
  // light
  spheres[5] = (Sphere) { .center = VEC3_xyz(0.0, 20, -30),
                          .radius = 3,
                          .radius2 = 3 * 3,
                          .surfaceColor = VEC3_xyz(0.00, 0.00, 0.00),
                          .emissionColor = VEC3_x(3),
                          .transparency = 0,
                          .reflection = 0 };
  render(spheres, num_spheres);

  return 0;
}
