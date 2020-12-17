#include "vec.h"

#include "assert.h"
#include "math.h"


/* --- VEC 2 --- */


vec2 vec2_new(float x, float y)
{
    vec2 v2 = (vec2){.x=x, .y=y};
    return v2;
}

vec2 vec2_add(const vec2 *a, const vec2 *b)
{
    return vec2_new(a->x + b->x, a->y + b->y);
}

vec2 vec2_sub(const vec2 *a, const vec2 *b)
{
    return vec2_new(a->x - b->x, a->y - b->y);
}

vec2 vec2_mul(const vec2 *v2, float by)
{
    return vec2_new(v2->x * by, v2->y * by);
}

vec2 vec2_div(const vec2 *v2, float by)
{
    assert(by != 0.0);
    return vec2_new(v2->x / by, v2->y / by);
}

float vec2_dot(const vec2 *a, const vec2 *b)
{
    return (a->x * b->x) + (a->y * b->y);
}

float vec2_len(const vec2 *v2)
{
    float dot = vec2_dot(v2, v2);
    assert(dot != 0.0f);
    return sqrtf(dot);
}

vec2 vec2_normal(const vec2 *v2)
{
    return vec2_div(v2, vec2_len(v2));
}

/* --- VEC 3 --- */

vec3 vec3_new(float x, float y, float z)
{
    vec3 v3 = (vec3){.x=x, .y=y, .z=z};
    return v3;
}

vec3 vec3_add(const vec3 *a, const vec3 *b)
{
    return vec3_new(a->x + b->x, a->y + b->y, a->z + b->z);
}

vec3 vec3_sub(const vec3 *a, const vec3 *b)
{
    return vec3_new(a->x - b->x, a->y - b->y, a->z - b->z);
}

vec3 vec3_mul(const vec3 *v3, float by)
{
    return vec3_new(v3->x * by, v3->y * by, v3->z * by);
}

vec3 vec3_div(const vec3 *v3, float by)
{
    assert(by != 0.0);
    return vec3_new(v3->x / by, v3->y / by, v3->z / by);
}

float vec3_dot(const vec3 *a, const vec3 *b)
{
    return (a->x * b->x) + (a->y * b->y) + (a->z  * b->z);
}

vec3 vec3_cross(const vec3 *a, const vec3 *b)
{
    float x = (a->y * b->z) - (a->z * b->y);
    float y = (a->z * b->x) - (a->x * b->z);
    float z = (a->x * b->y) - (a->y * b->x);
    return vec3_new(x, y, z);
}

float vec3_len(const vec3 *v3)
{
    float dot = vec3_dot(v3, v3);
    assert(dot != 0.0f);
    return sqrtf(dot);
}

vec3 vec3_normal(const vec3 *v3)
{
    return vec3_div(v3, vec3_len(v3));
}


/* --- VEC 4 --- */


vec4 vec4_new(float x, float y, float z, float a)
{
    vec4 v4 = (vec4){.x=x, .y=y, .z=z, .a=a};
    return v4;
}
