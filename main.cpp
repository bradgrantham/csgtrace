#include <array>
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

struct shading_info
{
    vec3f p;
    vec3f n;
    vec3f c;
    enum direction { IN, OUT } d;
};

struct shape
{
    virtual bool intersect(const ray& ray, const range& r, shading_info* shade) = 0;
    virtual vec3f random_from_point(const vec3f& point) = 0;
};

bool intersect_triangle(const vec3f* vertices, const ray& ray, const range& r, triangle_hit* hit)
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

struct triangle : public shape
{
    std::array<vec3f, 3> vertices;
    triangle(const vec3f* vertices_)
    {
        std::copy(vertices_, vertices_ + 3, vertices.begin());
    }
    virtual bool intersect(const ray& ray, const range& r, shading_info* shade)
    {
        triangle_hit hit;
        bool hits = intersect_triangle(vertices.data(), ray, r, &hit);
        if(hits) {
            shade->p = uvwTriangleToPoint(vertices.data(), hit.uvw[0], hit.uvw[1], hit.uvw[2]);
            shade->n = {0, 0, 0}; // XXX
            shade->c = {1, 1, 1}; // XXX
        }
        return hits;
    }
    virtual vec3f random_from_point(const vec3f& point)
    {
        // XXX ignore point for now
        return randomInTriangle(vertices.data());
    }
};

struct csg_node
{
    struct direction_hit
    {
        float t;
        shading_info::direction d;
        vec3f n;
    };
    virtual shading_info::direction classify(const vec3f& point) = 0;
    virtual bool intersect(const ray& ray, const range& r, direction_hit* hit) = 0;
    virtual vec3f random_from_point(const vec3f& from) = 0;
};

struct csg_sphere : csg_node
{
    vec3f center;
    float radius;

    csg_sphere(const vec3f &center, float radius) :
        center(center),
        radius(radius)
    {}

    virtual shading_info::direction classify(const vec3f& point)
    {
        return (vec_length(point - center) < radius) ? shading_info::IN : shading_info::OUT;
    }

    virtual bool intersect(const ray& ray, const range& r, direction_hit* hit)
    {
        vec3f oc = ray.m_origin - center;
        float a = vec_dot(ray.m_direction, ray.m_direction);
        float half_b = vec_dot(oc, ray.m_direction);
        float c = vec_dot(oc, oc) - radius * radius;
        float discriminant = half_b * half_b - a * c;
        if (discriminant < 0) {
            return false;
        }

        float sd = sqrtf(discriminant);

        // Find the nearest root that lies in the acceptable range.
        float root = (-half_b - sd) / a;
        if (root < r.t0 || r.t1 < root) {
            root = (-half_b + sd) / a;
            if (root < r.t0 || r.t1 < root) {
                return false;
            }
        }

        vec3f point = ray.at(root);
        vec3f normal = (point - center) / radius;

        hit->t = root;
        bool into = (vec_dot(ray.m_direction, normal) < 0);
        hit->d = into ? shading_info::IN : shading_info::OUT;
        hit->n = into ? normal : -normal;

        return true;
    }

    virtual vec3f random_from_point(const vec3f& from)
    {
        return {0, 0, 0}; // XXX XXX!
    }
};

#if 0

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

#endif

struct csg_shape : public shape
{
    csg_node* root;

    csg_shape(csg_node *root) :
        root(root)
    { }

    // mapping from coord to color, maybe texture, maybe shader

    virtual bool intersect(const ray& ray, const range& r, shading_info* shade)
    {
        csg_node::direction_hit hit;
        bool hits = root->intersect(ray, r, &hit);
        if(hits) {
            shade->p = ray.m_origin + ray.m_direction * hit.t;
            shade->n = hit.n;
            shade->c = {1, 1, 1}; // XXX
            shade->d = hit.d; // XXX
        }

        return hits;
    }

    virtual vec3f random_from_point(const vec3f& point)
    {
        // XXX ignore point for now
        return root->random_from_point(point);
    }
};


vec3f cast(float u, float v, shape* scene, shape* light)
{
    vec3f onplane = {(u - .5f) * 2, (v - .5f) * 2, -1.0f}; // one meter away down Z
    shading_info shade;

    ray blockerTestRay = {{0, 0, 0}, onplane};

    if(scene->intersect(blockerTestRay, {-FLT_MAX, FLT_MAX}, &shade)) {

        // return {1, 1, 1};
        vec3f lightpoint = light->random_from_point(onplane);
        return vec3f(1, 1, 1) * std::clamp(vec_dot(shade.n, vec_normalize(lightpoint - shade.p)), 0.0f, 1.0f);

    } else {

        /* cast to random points on light */
        vec3f lightpoint = light->random_from_point(onplane);

        ray shadowTestRay = {onplane, vec_normalize(lightpoint - onplane)};

        if(!scene->intersect(shadowTestRay, {-FLT_MAX, FLT_MAX}, &shade)) {

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
    constexpr int sampleCount = 16;
    uint8_t img[height * width * 3];

    vec3f light_vertices[3] = {{.3, .3, 0}, {.3, .6, 0}, {.6, .6, 0}};
    triangle light(light_vertices);

    vec3f blocker_vertices[3] = {{0, 0, -1}, {.4, .3, -.8}, {.3, .4, -.8}};
    triangle blocker(blocker_vertices);

    csg_sphere blocker_sphere_node { { .2, .2, -.8}, .2 };
    csg_shape blocker2(&blocker_sphere_node);

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

                sum += cast(u, v, &blocker2, &light);
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
