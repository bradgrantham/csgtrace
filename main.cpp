#include <cstdio>
#include <cassert>
#include <cstdint>
#include <cmath>
#include "vectormath.h"

struct range
{
    float t0;
    float t1;
};

struct triangle_hit
{
    float t;
    float uvw[3];
};

#if 0
struct shading_info
{
    vec3f p;
    vec3f n;
    vec3f c;
};

struct csg_node
{
    enum direction { IN, OUT } ;
    struct direction_hit
    {
        float t;
        direction what;
        vec3f n;
    };
    virtual direction classify(const vec3f& point) = 0;
    virtual bool intersect(const ray& ray, const range& r, direction_hit* hit) = 0;
};

struct csg_sphere : csg_node
{
};

struct csg_transformation : csg_node
{
};

// half_space
// cylinder
// cone
// torus
// polyhedron
// union
// intersection
// difference

struct csg_object
{
    struct surface_info
    {
        vec3f p;
        vec3f n;
    };
    csg_node* root;
    // mapping from coord to color, maybe texture, maybe shader, who knows
    // intersect returning a shading_info?
};

#endif

bool intersect(const vec3f* vertices, const ray& ray, const range& r, triangle_hit* hit)
{
    vec3f e0 = vertices[1] - vertices[0];
    vec3f e1 = vertices[2] - vertices[0];

    vec3f P = vec_cross(ray.m_direction, e1);

    float det = vec_dot(e0, P);

    const float epsilon = 0.0000001; // .000001 too small

    if(det > -epsilon && det < epsilon) {
        return false;
    }

    float inv_det = 1.0f / det;

    vec3f T = ray.m_origin - vertices[0];

    float u = vec_dot(T, P) * inv_det;
    if(u < 0.0f || u > 1.0f) {
        return false;
    }

    vec3f Q = vec_cross(T, e0);
    float v = vec_dot(ray.m_direction, Q) * inv_det;
    if(v < 0.0f || u + v > 1.0f) {
        return false;
    }

    float d = vec_dot(e1, Q) * inv_det;
    if(d < r.t0 || d > r.t1) {
        return false;
    }

    hit->t = d;
    hit->uvw[0] = 1 - u - v;
    hit->uvw[1] = u;
    hit->uvw[2] = v;

    return true;
}

vec3f uvwTriangleToPoint(const vec3f* vertices, float u, float v, float w)
{
    return vertices[0] * u + vertices[1] * v + vertices[2] * w;
}

vec3f randomInTriangle(const vec3f* vertices)
{
    float u = drand48();
    float v = (1 - u) * drand48();
    float w = 1 - u * v;
    return uvwTriangleToPoint(vertices, u, v, w);
}

vec3f cast(float u, float v)
{
    vec3f onplane = {(u - .5f) * 2, (v - .5f) * 2, -1.0f}; // one meter away down Z
    vec3f light[3] = {{.3, .3, 0}, {.3, .6, 0}, {.6, .6, 0}};
    vec3f blocker[3] = {{0, 0, -1}, {.4, .3, -.8}, {.3, .4, -.8}};

    triangle_hit hit;

    ray blockerTestRay = {{0, 0, 0}, onplane};

    if(intersect(blocker, blockerTestRay, {-FLT_MAX, FLT_MAX}, &hit)) {

        return {1, 1, 1};

    } else {

        /* cast to random points on light */
        vec3f lightpoint = randomInTriangle(light);

        ray shadowTestRay = {onplane, vec_normalize(lightpoint - onplane)};

        if(!intersect(blocker, shadowTestRay, {-FLT_MAX, FLT_MAX}, &hit)) {

            // return {u, v, .5};
            return {fabsf(shadowTestRay.m_direction[0]), fabsf(shadowTestRay.m_direction[1]), fabsf(shadowTestRay.m_direction[2])};
        }
    }

    return {0, 0, 0};
}

int main(int argc, char **argv)
{
    constexpr int width = 512;
    constexpr int height = 512;
    constexpr int sampleCount = 1000;
    uint8_t img[height * width * 3];

    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {

            uint8_t *pixel = img + (x + y * width) * 3;

            pixel[0] = 0.0f;
            pixel[1] = 0.0f;
            pixel[2] = 0.0f;

            vec3f sum = {0, 0, 0};

            for(int sample = 0; sample < sampleCount; sample++) {
                float u = (x + drand48()) / width;
                float v = (y + drand48()) / height;

                sum += cast(u, v);
            }

            pixel[0] = sum.x / sampleCount * 255;
            pixel[1] = sum.y / sampleCount * 255;
            pixel[2] = sum.z / sampleCount * 255;
        }
    }

    FILE *fp = fopen("output.ppm", "wb");
    assert(fp);
    fprintf(fp, "P6 %d %d 255\n", width, height);
    for(int y = 0; y < height; y++) {
        fwrite(img + 3 * width * (height - 1 - y), 3, width, fp);
    }
    fclose(fp);
}
