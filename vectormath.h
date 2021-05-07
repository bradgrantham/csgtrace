#ifndef __VECTORMATH_H__
#define __VECTORMATH_H__

#include <math.h>
// #include <values.h>
#include <float.h>

/*
Style to do:
    vec{234}{i,f} should be specialized from std::array?
    size_t, not int, in arrays?
    run through clang-format to fix single-line, no-brace code
*/

template <class V>
float vec_dot(const V& v0, const V& v1)
{
    float total = 0;
    for(int i = 0; i < V::dimension(); i++) total += v0[i] * v1[i];
    return total;
}

template <class V>
float vec_length(const V& v)
{
    float sum = 0;
    for(int i = 0; i < V::dimension(); i++) sum += v[i] * v[i];
    return sqrtf(sum);
}

template <class V>
float vec_length_sq(const V& v)
{
    float sum = 0;
    for(int i = 0; i < V::dimension(); i++) sum += v[i] * v[i];
    return sum;
}

template <class V>
V vec_normalize(const V& v)
{
    float l = vec_length(v);
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = v[i] / l;
    return tmp;
}

template <class V>
V vec_blend(const V& v0, float w0, const V& v1, float w1)
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = v0[i] * w0 + v1[i] * w1;
    return tmp;
}

template <class V>
V vec_scale(const V& v0, float w0)
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = v0[i] * w0;
    return tmp;
}

/* normalized i, n */
/* doesn't normalize r */
/* r = u - 2 * n * dot(n, u) */
template <class V>
V vec_reflect(const V& i, const V& n)
{
    V tmp;
    tmp = vec_blend(i, 1.0f, n, -2.0f * vec_dot(i, n));
    return tmp;
}

template <class V>
bool vec_refract(float eta, const V& i, const V& n, V& result)
{
    float dot = vec_dot(n, i); 

    float k = 1.0f - eta * eta * (1.0f - dot * dot);

    if(k < 0.0f) {
        result = V(0, 0, 0);
        return true;
    }

    result = eta * i - (eta * dot + sqrtf(k)) * n;
    return false;
}

#if 0
template <class V>
bool operator==(const V& v0, const V& v1)
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++)
        if(v0[i] != v1[i])
            return false;
    return true;
}
#endif
    
template <class V>
V operator+(const V& v0, const V& v1)
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = v0[i] + v1[i];
    return tmp;
}

template <class V>
V operator-(const V& v0, const V& v1)
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = v0[i] - v1[i];
    return tmp;
}

template <class V>
V operator-(const V& v)
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = -v[i];
    return tmp;
}

template <class V>
V operator*(float w, const V& v) 
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = v[i] * w;
    return tmp;
}

template <class V>
V operator/(float w, const V& v)
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = v[i] / w;
    return tmp;
}

template <class V>
V operator*(const V& v, float w)
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = v[i] * w;
    return tmp;
}

template <class V>
V operator/(const V& v, float w)
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = v[i] / w;
    return tmp;
}

template <class V>
V operator*(const V& v0, const V& v1)
{
    V tmp;
    for(int i = 0; i < V::dimension(); i++) tmp[i] = v0[i] * v1[i];
    return tmp;
}

struct vec2f
{
    float m_v[2];
    float &x, &y;
    static constexpr int dimension() { return 2; }
    typedef float comp_type;

    vec2f(void) :
        x(m_v[0]),
        y(m_v[1])
    { }

    vec2f(float x, float y) :
        x(m_v[0]),
        y(m_v[1])
    { set(x, y); }

    void set(float x, float y)
        { m_v[0] = x; m_v[1] = y; }

    vec2f(float v) :
        x(m_v[0]),
        y(m_v[1])
    {
	for(int i = 0; i < 2; i++) m_v[i] = v;
    }

    vec2f(const float *v) :
        x(m_v[0]),
        y(m_v[1])
    {
	for(int i = 0; i < 2; i++) m_v[i] = v[i];
    }

    vec2f(const vec2f &v) :
        x(m_v[0]),
        y(m_v[1])
    {
	for(int i = 0; i < 2; i++) m_v[i] = v[i];
    }

    vec2f &operator=(const float *v) {
	for(int i = 0; i < 2; i++) m_v[i] = v[i];
	return *this;
    }

    vec2f &operator=(const vec2f& v) {
	for(int i = 0; i < 2; i++) m_v[i] = v[i];
	return *this;
    }

    vec2f(float *v) :
        x(m_v[0]),
        y(m_v[1])
        { set(v); }

    operator const float*() const { return m_v; }
    operator float*() { return m_v; }

    float& operator[] (int i)
        { return m_v[i]; }

    const float& operator[] (int i) const
        { return m_v[i]; }

    void clear() { 
	for(int i = 0; i < 2; i++) m_v[i] = 0;
    }

    void set(const float *v)
	{ for(int i = 0; i < 2; i++) m_v[i] = v[i]; }

    vec2f& normalize() {
	*this = vec_normalize(*this);
	return *this;
    }

    vec2f operator*=(float w) {
	for(int i = 0; i < 2; i++) m_v[i] *= w;
	return *this;
    }

    vec2f operator/=(float w) {
	for(int i = 0; i < 2; i++) m_v[i] /= w;
	return *this;
    }

    vec2f operator+=(const vec2f& v) {
	for(int i = 0; i < 2; i++) m_v[i] += v[i];
	return *this;
    }

    vec2f operator-=(const vec2f& v) {
	for(int i = 0; i < 2; i++) m_v[i] -= v[i];
	return *this;
    }
};

struct vec3f
{
    float m_v[3];
    float &x, &y, &z;
    static constexpr int dimension() { return 3; }
    typedef float comp_type;

    vec3f(void) :
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2])
    { }

    void set(float x, float y, float z)
        { m_v[0] = x; m_v[1] = y; m_v[2] = z;}

    vec3f(float x, float y, float z) :
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2])
        { set(x, y, z); }

    vec3f cross(const vec3f& v1) {
	vec3f tmp;
	tmp[0] = m_v[1] * v1[2] - m_v[2] * v1[1];
	tmp[1] = m_v[2] * v1[0] - m_v[0] * v1[2];
	tmp[2] = m_v[0] * v1[1] - m_v[1] * v1[0];
	*this = tmp;
	return *this;
    }

    vec3f(float v) :
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2])
    {
	for(int i = 0; i < 3; i++) m_v[i] = v;
    }

    vec3f(const float *v) :
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2])
    {
	for(int i = 0; i < 3; i++) m_v[i] = v[i];
    }

    vec3f(const vec3f &v) :  
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2])
    {
	for(int i = 0; i < 3; i++) m_v[i] = v[i];
    }

    vec3f &operator=(const vec3f& v) {
	for(int i = 0; i < 3; i++) m_v[i] = v[i];
	return *this;
    }

    vec3f &operator=(float v) {
	for(int i = 0; i < 3; i++) m_v[i] = v;
	return *this;
    }

    vec3f(float *v) : 
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2])
    { set(v); }

    operator const float*() const { return m_v; }
    operator float*() { return m_v; }

    float& operator[] (int i)
        { return m_v[i]; }

    const float& operator[] (int i) const
        { return m_v[i]; }

    void clear() { 
	for(int i = 0; i < 3; i++) m_v[i] = 0;
    }

    void set(const float *v)
	{ for(int i = 0; i < 3; i++) m_v[i] = v[i]; }

    float length() const {
	float sum = 0;
	for(int i = 0; i < 3; i++) sum += m_v[i] * m_v[i];
	return (float)sqrtf(sum);
    }

    vec3f& normalize() {
	*this = vec_normalize(*this);
	return *this;
    }

    vec3f operator*=(float w) {
	for(int i = 0; i < 3; i++) m_v[i] *= w;
	return *this;
    }

    vec3f operator/=(float w) {
	for(int i = 0; i < 3; i++) m_v[i] /= w;
	return *this;
    }

    vec3f operator+=(const vec3f& v) {
	for(int i = 0; i < 3; i++) m_v[i] += v[i];
	return *this;
    }

    vec3f operator-=(const vec3f& v) {
	for(int i = 0; i < 3; i++) m_v[i] -= v[i];
	return *this;
    }

    vec3f operator*=(const vec3f& v) {
	for(int i = 0; i < 3; i++) m_v[i] *= v[i];
	return *this;
    }

    vec3f operator/=(const vec3f& v) {
	for(int i = 0; i < 3; i++) m_v[i] /= v[i];
	return *this;
    }
};

vec3f vec_cross(const vec3f& v0, const vec3f& v1)
{
    vec3f tmp;
    tmp[0] = v0[1] * v1[2] - v0[2] * v1[1];
    tmp[1] = v0[2] * v1[0] - v0[0] * v1[2];
    tmp[2] = v0[0] * v1[1] - v0[1] * v1[0];
    return tmp;
}

struct vec4f
{
    float m_v[4];
    float &x, &y, &z, &w;
    static constexpr int dimension() { return 4; }
    typedef float comp_type;

    vec4f(void) :
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2]),
        w(m_v[3])
    { }

    void set(float x, float y, float z, float w)
    {
        m_v[0] = x; m_v[1] = y; m_v[2] = z; m_v[3] = w;
    }

    vec4f(float x, float y, float z, float w) :
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2]),
        w(m_v[3])
    { set(x, y, z, w); }

    vec4f(float v)  :
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2]),
        w(m_v[3])
    {
	for(int i = 0; i < 4; i++) m_v[i] = v;
    }

    vec4f(const float *v) :
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2]),
        w(m_v[3])
    {
	for(int i = 0; i < 4; i++) m_v[i] = v[i];
    }

    vec4f(const vec4f &v) :
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2]),
        w(m_v[3])
    {
	for(int i = 0; i < 4; i++) m_v[i] = v[i];
    }

    vec4f &operator=(const vec4f& v) {
	for(int i = 0; i < 4; i++) m_v[i] = v[i];
	return *this;
    }

    vec4f(float *v) :
        x(m_v[0]),
        y(m_v[1]),
        z(m_v[2]),
        w(m_v[3])
    {
        set(v);
        }

    operator const float*() const { return m_v; }
    operator float*() { return m_v; }

    float& operator[] (int i)
        { return m_v[i]; }

    const float& operator[] (int i) const
        { return m_v[i]; }

    void clear() { 
	for(int i = 0; i < 4; i++) m_v[i] = 0;
    }

    void set(const float *v)
	{ for(int i = 0; i < 4; i++) m_v[i] = v[i]; }

    float length() const {
	float sum = 0;
	for(int i = 0; i < 4; i++) sum += m_v[i] * m_v[i];
	return (float)sqrtf(sum);
    }

    vec4f& normalize() {
	*this = vec_normalize(*this);
	return *this;
    }

    vec4f operator*=(float w) {
	for(int i = 0; i < 4; i++) m_v[i] *= w;
	return *this;
    }

    vec4f operator/=(float w) {
	for(int i = 0; i < 4; i++) m_v[i] /= w;
	return *this;
    }

    vec4f operator+=(const vec4f& v) {
	for(int i = 0; i < 4; i++) m_v[i] += v[i];
	return *this;
    }

    vec4f operator-=(const vec4f& v) {
	for(int i = 0; i < 4; i++) m_v[i] -= v[i];
	return *this;
    }
};

//
// With your left hand up, fingers up, palm facing away, thumb facing to
// the right, thumb is v0-v1, index finger is v0-v2 : plane normal
// sticks out the back of your the hand towards you.
//
inline vec4f make_plane(const vec3f& v0, const vec3f& v1, const vec3f& v2)
{
    vec3f xaxis, yaxis;
    vec3f plane;

    xaxis = vec_blend(v1, 1.0, v0, -1.0);

    yaxis = vec_blend(v2, 1.0, v0, -1.0);

    plane = vec_cross(xaxis, yaxis);
    plane.normalize();

    float D = vec_dot(-v0, plane);
    return vec4f(plane[0], plane[1], plane[2], D);
}

struct rot4f : public vec4f
{
    void set_axis(float x, float y, float z) {
	m_v[1] = x;
	m_v[2] = y;
	m_v[3] = z;
    }

    void set_axis(vec3f &axis) {
	m_v[1] = axis[0];
	m_v[2] = axis[1];
	m_v[3] = axis[2];
    }

    rot4f& mult(const rot4f& m1, const rot4f &m2);
};

rot4f operator*(const rot4f& r1, const rot4f& r2);

struct mat4f
{
    float m_v[16];
    static constexpr int dimension() { return 16; }
    typedef float comp_type;
    static mat4f identity;

    mat4f() { }

    mat4f(float m00, float m01, float m02, float m03,
	float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23,
	float m30, float m31, float m32, float m33) {

	m_v[0] = m00; m_v[1] = m01; m_v[2] = m02; m_v[3] = m03;
	m_v[4] = m10; m_v[5] = m11; m_v[6] = m12; m_v[7] = m13;
	m_v[8] = m20; m_v[9] = m21; m_v[10] = m22; m_v[11] = m23;
	m_v[12] = m30; m_v[13] = m31; m_v[14] = m32; m_v[15] = m33;
    }

    // mat4f::mult_nm does not perform inverse transpose - just multiplies with
    //     v[3] = 0
    vec3f mult_nm(vec3f &in) {
	int i;
	vec4f t;

	for(i = 0; i < 4; i++)
	    t[i] =
		m_v[0 + i] * in[0] + 
		m_v[4 + i] * in[1] + 
		m_v[8 + i] * in[2];

	t[0] /= t[3];
	t[1] /= t[3];
	t[2] /= t[3];
	return vec3f(t[0], t[1], t[2]);
    }

    mat4f& transpose(mat4f& in) {
	mat4f t;
	int i, j;

	t = in;
	for(i = 0; i < 4; i++)
	    for(j = 0; j < 4; j++) 
		m_v[i + j * 4] = t[j + i * 4];

	return *this;
    }

    float determinant() const {
	return (m_v[0] * m_v[5] - m_v[1] * m_v[4]) *
	    (m_v[10] * m_v[15] - m_v[11] * m_v[14]) + 
	    (m_v[2] * m_v[4] - m_v[0] * m_v[6]) *
	    (m_v[9] * m_v[15] - m_v[11] * m_v[13]) + 
	    (m_v[0] * m_v[7] - m_v[3] * m_v[4]) *
	    (m_v[9] * m_v[14] - m_v[10] * m_v[13]) + 
	    (m_v[1] * m_v[6] - m_v[2] * m_v[5]) *
	    (m_v[8] * m_v[15] - m_v[11] * m_v[12]) + 
	    (m_v[3] * m_v[5] - m_v[1] * m_v[7]) *
	    (m_v[8] * m_v[14] - m_v[10] * m_v[12]) + 
	    (m_v[2] * m_v[7] - m_v[3] * m_v[6]) *
	    (m_v[8] * m_v[13] - m_v[9] * m_v[12]);
    }

    bool invert(const mat4f& in, bool singular_fail = true);
    bool invert() { return invert(*this); }

    static mat4f translation(float x, float y, float z) {
	mat4f m(identity);
	m[12] = x;
	m[13] = y;
	m[14] = z;

	return m;
    }

    static mat4f scale(float x, float y, float z) {
	mat4f m(identity);
	m[0] = x;
	m[5] = y;
	m[10] = z;

	return m;
    }

    static mat4f rotation(float a, float x, float y, float z) {
	mat4f m;
	float c, s, t;

	c = (float)cos(a);
	s = (float)sin(a);
	t = 1.0f - c;

	m[0] = t * x * x + c;
	m[1] = t * x * y + s * z;
	m[2] = t * x * z - s * y;
	m[3] = 0;

	m[4] = t * x * y - s * z;
	m[5] = t * y * y + c;
	m[6] = t * y * z + s * x;
	m[7] = 0;

	m[8] = t * x * z + s * y;
	m[9] = t * y * z - s * x;
	m[10] = t * z * z + c;
	m[11] = 0;

	m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;

	return m;
    }

    mat4f(const rot4f& r) {
	(*this) = rotation(r[0], r[1], r[2], r[3]);
    }

    void calc_rot4f(rot4f *out) const;

    mat4f& mult(mat4f& m1, mat4f &m2) {
	mat4f t;
	int i, j;

	for(j = 0; j < 4; j++)
	    for(i = 0; i < 4; i++)
	       t[i * 4 + j] = m1[i * 4 + 0] * m2[0 * 4 + j] +
		   m1[i * 4 + 1] * m2[1 * 4 + j] +
		   m1[i * 4 + 2] * m2[2 * 4 + j] +
		   m1[i * 4 + 3] * m2[3 * 4 + j];

	*this = t;
	return *this;
    }

    mat4f(const float *v) {
	for(int i = 0; i < 16; i++) m_v[i] = v[i];
    }

    mat4f(const mat4f &v) {
	for(int i = 0; i < 16; i++) m_v[i] = v[i];
    }

    mat4f &operator=(const mat4f& v) {
	for(int i = 0; i < 16; i++) m_v[i] = v[i];
	return *this;
    }

    mat4f(float *v)
        { set(v); }

    float& operator[] (int i)
        { return m_v[i]; }

    const float& operator[] (int i) const
        { return m_v[i]; }

    void clear() { 
	for(int i = 0; i < 16; i++) m_v[i] = 0;
    }

    void set(const float *v)
	{ for(int i = 0; i < 16; i++) m_v[i] = v[i]; }

    void store(float *v)
	{ for(int i = 0; i < 16; i++) v[i] = m_v[i]; }
};

inline mat4f operator*(const mat4f& m1, const mat4f& m2)
{
    mat4f t;
    int i, j;

    for(j = 0; j < 4; j++)
	for(i = 0; i < 4; i++)
	   t[i * 4 + j] =
	       m1[i * 4 + 0] * m2[0 * 4 + j] +
	       m1[i * 4 + 1] * m2[1 * 4 + j] +
	       m1[i * 4 + 2] * m2[2 * 4 + j] +
	       m1[i * 4 + 3] * m2[3 * 4 + j];

    return t;
}

inline vec4f operator*(const vec4f& in, const mat4f& m)
{
    int i;
    vec4f t;

    for(i = 0; i < 4; i++)
	t[i] =
	    m[0 + i] * in[0] + 
	    m[4 + i] * in[1] + 
	    m[8 + i] * in[2] + 
	    m[12 + i] * in[3];
    return t;
}

#if 0
inline vec3f operator*(const vec3f& in, const mat4f& m)
{
    int i;
    vec3f t;

    for(i = 0; i < 3; i++)
	t[i] =
	    m[0 + i] * in[0] + 
	    m[4 + i] * in[1] + 
	    m[8 + i] * in[2] + 
	    m[12 + i];
    return t;
}
#else
inline vec3f operator*(const vec3f& in, const mat4f& m)
{
    int i;
    vec4f t;

    for(i = 0; i < 4; i++)
	t[i] =
	    m[0 + i] * in[0] + 
	    m[4 + i] * in[1] + 
	    m[8 + i] * in[2] + 
	    m[12 + i];
    return vec3f(t.m_v);
}
#endif

enum axis_t { X_AXIS = 0, Y_AXIS = 1, Z_AXIS = 2};

struct segment
{
    vec3f m_v0;
    vec3f m_v1;
    segment(const vec3f &v0, const vec3f &v1) :
        m_v0(v0),
        m_v1(v1)
    {}
};

struct ray
{
    vec3f m_origin;
    vec3f m_direction;

    ray(const vec3f& o, const vec3f& d) : m_origin(o), m_direction(d)
    { }

    ray(const segment &s) :
        m_origin(s.m_v0),
        m_direction(s.m_v1 - s.m_v0)
    { }

    ray() {}

    float length() const
    {
        return m_direction.length();
    }

    float at(int axis, float plane) const
    {
	if(m_direction[axis] > -.00001f && m_direction[axis] < 0.0f)
	    return -FLT_MAX;
	if(m_direction[axis] >= 0.0f && m_direction[axis] < 0.00001f)
	    return FLT_MAX;
	return (plane - m_origin[axis]) / m_direction[axis];
    }

    vec3f at(float t) const
    {
        return m_origin + m_direction * t;
    }
};

// Probably should do this by passing in inverse-transpose for direction
inline ray operator*(const ray& r, const mat4f& m)
{
    vec3f newo = r.m_origin * m;
    vec3f newd = (r.m_direction + r.m_origin) * m - newo;
    return ray(newo, newd);
}

inline float operator*(const ray& r, const vec4f& plane)
{
    vec3f normal(plane[0], plane[1], plane[2]);

    float factor = vec_dot(r.m_direction, normal);
    return (plane[3] - vec_dot(normal, r.m_origin)) / factor;
}

inline void transform_ray(const vec3f& origin, const vec3f& direction, const mat4f& m, vec3f *neworigin, vec3f *newdirection)
{
    vec3f oldorigin = origin; // in case in and out are the same

    *neworigin = oldorigin * m;
    *newdirection = (direction + oldorigin) * m - *neworigin;
}

inline void transform_ray(const vec4f& origin, const vec4f& direction, const mat4f& m, vec4f *neworigin, vec4f *newdirection)
{
    vec4f oldorigin = origin; // in case in and out are the same

    *neworigin = oldorigin * m;
    *newdirection = (oldorigin + direction) * m - *neworigin;
}

inline void transform_ray(vec3f* origin, vec3f* direction, const mat4f& m)
{
    vec3f neworigin, newdirection;

    neworigin = *origin * m;
    newdirection = (*direction + *origin) * m - neworigin;

    *origin = neworigin;
    *direction = newdirection;
}

#endif /* __VECTORMATH_H__ */
